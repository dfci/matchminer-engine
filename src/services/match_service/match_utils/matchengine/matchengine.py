import json
import logging
import networkx as nx

from src.utilities import settings as s
from src.data_store import key_names as kn
from src.data_store.validator import SamplesValidator
from src.services.match_service.match_utils.sort import Sort
from src.utilities.utilities import get_db, format_match_tree_code
from src.data_store.trial_matches_data_model import trial_matches_schema
from src.services.match_service.match_utils.matchengine.assess_node_utils import AssessNodeUtils
from src.services.match_service.match_utils.matchengine.intersect_results_utils import IntersectResultsUtils

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s', )


class MatchEngine(AssessNodeUtils, IntersectResultsUtils):

    def __init__(self, match_tree, trial_info):
        AssessNodeUtils.__init__(self)
        IntersectResultsUtils.__init__(self)

        self.trial_info = trial_info
        self.match_tree = match_tree
        self.diagnosis_level = None
        self.match_tree_nx = None
        self.sample_ids = None
        self.db = get_db()

        self.matches = None
        self.trial_matches = None
        self.validator = SamplesValidator(trial_matches_schema)

    def convert_match_tree_to_digraph(self):
        """
        Converts a match tree in JSON format into a networkx directed graph
        Given json object of MATCH clause, the function returns a directed graph

        :return: {diGraph object}
        """
        key = self.match_tree.keys()[0]
        value = self.match_tree[key]
        global_node = 1
        self.match_tree_nx = nx.DiGraph()
        s = []
        s.append([0, global_node, key, value])
        while len(s) > 0:
            current = s.pop(0)
            parent = current[0]
            node = current[1]
            key = current[2]
            value = current[3]
            self.match_tree_nx.add_node(node)
            self.match_tree_nx.add_edges_from([(parent, node)])
            self.match_tree_nx.node[node]['type'] = key
            if isinstance(value, dict):
                self.match_tree_nx.node[node]['value'] = value
            elif isinstance(value, list):
                for i in range(0, len(value)):
                    global_node += 1
                    s.append([node, global_node, value[i].keys()[0], value[i][value[i].keys()[0]]])

        self.match_tree_nx.remove_node(0)

    def traverse_match_tree(self):
        """
        Construct a MongoDB query that represents the match tree

        :return: {null}
        """
        for node_id in list(nx.dfs_postorder_nodes(self.match_tree_nx, source=1)):

            # access node and its children
            node = self.match_tree_nx.node[node_id]
            children = [self.match_tree_nx.node[n] for n in self.match_tree_nx.successors(node_id)]

            # clinical nodes
            if node['type'] == 'clinical':
                node = self.assess_clinical_node(node=node)
                node['matches'] = self._search_for_matching_records(node=node, sample_ids=self.sample_ids)

            # genomic nodes
            elif node['type'] == 'genomic':
                node = self.assess_genomic_node(node=node)
                node['matches'] = self._search_for_matching_records(node=node, sample_ids=self.sample_ids)

            # join child queries with "and"
            elif node['type'] == 'and' or node['type'] == 'or':
                self.intersect_results(node=node, children=children)

            else:
                raise ValueError('match tree node must be of type "clinical", "genomic", "and", or "or')

        self.matches = self.match_tree_nx.node[1]['matches']

    def _search_for_matching_records(self, node, sample_ids=None):
        """
        Search for any sample records that match the constructed query

        :param node: {digraph node}
        :param sample_ids: {list of str or null}
        :return: {null}
        """
        # include inclusion reasons in projection
        if 'genomic_inclusion_reasons' in node:
            proj = node['genomic_inclusion_reasons']
        elif 'clinical_inclusion_reasons' in node:
            proj = node['clinical_inclusion_reasons']
        else:
            proj = self.proj.copy()

        # subset patient data if applicable
        if sample_ids is not None:
            node['query'][kn.sample_id_col] = {'$in': sample_ids}

        # perform query
        matches = list(self.db[s.sample_collection_name].find(node['query'], proj))

        # add exclusion reasons to match results
        for match in matches:
            if kn.mutation_list_col in match:
                for variant in match[kn.mutation_list_col]:
                    if 'variant_level' in node:
                        variant[kn.mr_reason_level_col] = node['variant_level']

            elif kn.genomic_exclusion_reasons_col in node:
                match[kn.genomic_exclusion_reasons_col] = [node[kn.genomic_exclusion_reasons_col]]
            elif kn.clinical_exclusion_reasons_col in node:
                match[kn.clinical_exclusion_reasons_col] = [node[kn.clinical_exclusion_reasons_col]]

        return matches

    def create_trial_match_records(self):
        """
        Create trial match records from matched samples,
        including clinical and genomic reasons for each match

        :return: {null}
        """
        for sample in self.matches:

            sample[kn.tm_sort_order_col] = 0
            sample[kn.tm_trial_protocol_no_col] = self.trial_info['protocol_no']
            sample[kn.tm_trial_accrual_status_col] = self.trial_info['accrual_status']
            sample[kn.mr_trial_level_col] = self.trial_info['level']
            sample[kn.mr_trial_step_code_col] = self.trial_info['step_code']
            sample[kn.mr_trial_arm_code_col] = self.trial_info['arm_code'] if 'arm_code' in self.trial_info else None
            sample[kn.mr_trial_dose_code_col] = self.trial_info['dose_code'] if 'dose_code' in self.trial_info else None
            sample[kn.mr_coordinating_center_col] = self.trial_info[s.trial_coordinating_center_col]

            if not self.validator.validate_document(sample):
                raise ValueError('%s sample did not pass data validation: %s' % (sample[kn.sample_id_col],
                                                                                 self.validator.errors))

        # todo sort
        # trial_matches_df = self.sort_trial_matches()

        # todo add versioning
        if len(self.matches) > 0:
            query = {kn.tm_trial_protocol_no_col == self.trial_info['protocol_no']}
            # records = json.loads(trial_matches_df.T.to_json()).values()
            self.db[s.trial_match_collection_name].remove(query)
            res = self.db[s.trial_match_collection_name].insert_many(self.matches)
            logging.info('%s | %s | %d trial matches added' % (
                self.trial_info['protocol_no'],
                format_match_tree_code(step_code=self.trial_info['step_code'],
                                       arm_code=self.trial_info['arm_code'],
                                       dose_code=self.trial_info['dose_code']),
                len(res.inserted_ids)))
        else:
            logging.info('%s | %s | No trial matches' % (
                self.trial_info['protocol_no'],
                format_match_tree_code(step_code=self.trial_info['step_code'],
                                       arm_code=self.trial_info['arm_code'],
                                       dose_code=self.trial_info['dose_code'])))

    def sort_trial_matches(self):
        """
        Sort trial matches.

        :return: {Pandas dataframe}
        """
        if self.matches is None:
            return

        logging.info('Sorting trial matches')

        print '---debug---'
        print self.matches[0]
        print
        sort = Sort(trial_matches=self.matches)
        return sort.add_sort_order()
