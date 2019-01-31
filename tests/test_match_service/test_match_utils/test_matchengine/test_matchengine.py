from src.utilities import settings as s
from src.data_store import key_names as kn
from tests.test_match_service import TestQueryUtilitiesShared
from src.services.match_service.match_utils.trial_utils import TrialUtils
from src.services.match_service.match_utils.matchengine.matchengine import MatchEngine
from src.data_store.trial_matches_data_model import trial_matches_schema

import datetime as dt


class TestMatchEngine(TestQueryUtilitiesShared):

    def setUp(self):
        super(TestMatchEngine, self).setUp()

        self.db.testSamples.drop()
        self.db[s.sample_collection_name].drop()
        self.db[s.trial_match_collection_name].drop()
        self.db.testSamples.insert_many(self.test_cases)

    def tearDown(self):
        self.db.testSamples.drop()
        self.db[s.sample_collection_name].drop()
        self.db[s.trial_match_collection_name].drop()

    @staticmethod
    def _set_up_matchengine(trial):
        t = TrialUtils(trial=trial)
        return t.parse_match_trees_from_trial()[0]

    def test_complex_conversion(self):

        # set up matching records
        doc1 = {
            kn.sample_id_col: 'DEV-01',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'EGFR', kn.transcript_exon_col: 19, kn.variant_class_col: 'In_Frame_Del'}
            ]
        }  # (y)
        doc2 = {
            kn.sample_id_col: 'DEV-02',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.XXXX'
                },
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                }
            ]
        }  # (y)
        doc3 = {
            kn.sample_id_col: 'DEV-03',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.XXXX'
                },
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                }
            ],
            kn.sv_list_col: [{kn.sv_comment_col: 'This contains an ALK SV'}]
        }  # (n)
        doc4 = {
            kn.sample_id_col: 'DEV-04',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Melanoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.XXXX'
                },
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                }
            ]
        }  # (n)
        doc5 = {
            kn.sample_id_col: 'DEV-05',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.L833V'
                },
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                },
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.V600D'
                }
            ]
        }  # (y)
        doc6 = {
            kn.sample_id_col: 'DEV-06',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.XXXX'
                },
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                },
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.V600E'
                }
            ]
        }  # (y)

        self.db[s.sample_collection_name].insert_many([doc1, doc2, doc3, doc4, doc5, doc6])

        # find matches
        me = MatchEngine(match_tree=self.complex_match_tree, trial_info={}, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.convert_match_tree_to_digraph()
        me.traverse_match_tree()
        self._print(me.matches)
        # todo is DEV-05 missing the exclusion criteria?
        assert sorted(list(set([i[kn.sample_id_col] for i in me.matches]))) == ['DEV-01', 'DEV-02', 'DEV-05'], \
            sorted([i[kn.sample_id_col] for i in me.matches])

    def test_combination_match(self):

        # set up matching records
        doc1 = {
            kn.sample_id_col: 'DEV-01',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                },
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.V600E'
                }
            ]
        }
        doc2 = {
            kn.sample_id_col: 'DEV-02',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.V600E'
                }
            ]
        }
        doc3 = {
            kn.sample_id_col: 'DEV-03',
            kn.mrn_col: '01',
            kn.vital_status_col: 'alive',
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'EGFR',
                    kn.variant_category_col: s.variant_category_mutation_val,
                    kn.protein_change_col: 'p.E709K'
                }
            ]
        }

        self.db[s.sample_collection_name].insert_many([doc1, doc2, doc3])

        # find matches
        me = MatchEngine(match_tree=self.simple_and_tree, trial_info={}, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.convert_match_tree_to_digraph()
        me.traverse_match_tree()
        print '---AND---'
        self._print(me.matches)

        me = MatchEngine(match_tree=self.simple_or_tree, trial_info={}, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.convert_match_tree_to_digraph()
        me.traverse_match_tree()
        print '\n---OR---'
        self._print(me.matches)

        me = MatchEngine(match_tree=self.complex_or_tree, trial_info={}, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.convert_match_tree_to_digraph()
        me.traverse_match_tree()
        print '\n---COMPLEX OR---'
        self._print(me.matches)

    def test_convert_match_tree_to_digraph(self):

        me = MatchEngine(match_tree=self.simple_mutation_match_tree, trial_info={}, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.convert_match_tree_to_digraph()
        assert me.match_tree_nx is not None
        assert me.match_tree_nx.node[1] == {'type': 'and'}
        assert me.match_tree_nx.node[2] == {'type': 'genomic', 'value': self.simple_mutation_match_tree['and'][0]['genomic']}
        assert me.match_tree_nx.node[3] == {'type': 'clinical', 'value': self.simple_mutation_match_tree['and'][1]['clinical']}

    def test_search_for_matching_records(self):

        self.db[s.sample_collection_name].insert_many([
            self.test_case_lung,
            self.test_case_colon,
            self.test_case_braf_v600e
        ])

        # clinical inclusion
        node = {
            'type': 'clinical',
            'query': {kn.oncotree_primary_diagnosis_name_col: 'Lung'},
            'clinical_inclusion_reasons': {
                '_id': 0, kn.sample_id_col: 1,
                kn.oncotree_primary_diagnosis_name_col: 1
            }
        }
        me = MatchEngine(match_tree=None, trial_info=None, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        m = me._search_for_matching_records(node=node)
        assert m == [{
            kn.sample_id_col: 'TEST-SAMPLE-LUNG',
            kn.oncotree_primary_diagnosis_name_col: 'Lung'
        }], m

        # clinical exclusion
        node = {
            'type': 'clinical',
            'query': {kn.oncotree_primary_diagnosis_name_col: {'$ne': 'Lung'}},
            kn.clinical_exclusion_reasons_col: {kn.oncotree_primary_diagnosis_name_col: 'Lung'}
        }
        me = MatchEngine(match_tree=None, trial_info=None, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        m = me._search_for_matching_records(node=node)
        assert m[0] == {
            kn.sample_id_col: 'TEST-SAMPLE-COLON',
            kn.clinical_exclusion_reasons_col: [{
                kn.oncotree_primary_diagnosis_name_col: 'Lung'
            }]
        }, m[0]

        # genomic inclusion
        node = {
            'type': 'genomic',
            'query': {kn.mutation_list_col: {'$elemMatch': {kn.hugo_symbol_col: 'BRAF'}}},
            'variant_level': 'gene',
            'genomic_inclusion_reasons': {
                '_id': 0, kn.sample_id_col: 1,
                kn.oncotree_primary_diagnosis_name_col: 1,
                '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 1,
                '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 1,
                '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 1,
                '%s.%s' % (kn.mutation_list_col, kn.cnv_call_col): 1,
                '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 1,
            },
            'unwind': '$%s' % kn.mutation_list_col,
            'match_reason': {'%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF'},
            'include': True
        }
        me = MatchEngine(match_tree=None, trial_info=None, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        m = me._search_for_matching_records(node=node)
        assert m == [{
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-V600E',
            kn.oncotree_primary_diagnosis_name_col: 'Breast',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.ref_residue_col: 'p.V600',
                kn.mr_reason_level_col: 'gene',
            }]
        }], m

        # genomic exclusion
        node = {
            'type': 'genomic',
            'query': {kn.mutation_list_col: {'$not': {'$elemMatch': {kn.hugo_symbol_col: 'BRAF'}}}},
            'variant_level': 'gene',
            kn.genomic_exclusion_reasons_col: {
                kn.hugo_symbol_col: 'BRAF',
                kn.mutation_list_col: s.variant_category_mutation_val
            }
        }
        me = MatchEngine(match_tree=None, trial_info=None, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        m = me._search_for_matching_records(node=node)
        assert m[1] == {
            kn.sample_id_col: 'TEST-SAMPLE-COLON',
            kn.genomic_exclusion_reasons_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.mutation_list_col: s.variant_category_mutation_val
            }]
        }, m[1]

    def test_create_trial_match_records(self):

        trial_info = {
            'protocol_no': '00-000',
            'accrual_status': 'open',
            'level': 'step',
            'step_code': 'STEP CODE',
            'arm_code': None
        }
        me = MatchEngine(match_tree=None, trial_info=trial_info, mongo_uri='mongodb://localhost:27017', mongo_dbname='matchminer')
        me.matches = [
            {
                kn.genomic_exclusion_reasons_col: [
                    {
                        kn.variant_category_col: s.variant_category_sv_val,
                        kn.hugo_symbol_col: 'ALK'
                    },
                ],
                kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
                kn.mr_diagnosis_level_col: '_solid_',
                kn.mutation_list_col: [{kn.hugo_symbol_col: 'BRAF', 'level': 'gene'}],
                kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
                kn.sample_id_col: 'DEV-01',
                kn.mrn_col: '01',
                kn.vital_status_col: 'alive'
            },
            {
                kn.genomic_exclusion_reasons_col: [
                    {
                        kn.variant_category_col: s.variant_category_sv_val,
                        kn.hugo_symbol_col: 'ALK'
                    },
                ],
                kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
                kn.mr_diagnosis_level_col: 'specific',
                kn.mutation_list_col: [{kn.hugo_symbol_col: 'EGFR', 'level': 'gene'}],
                kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
                kn.sample_id_col: 'DEV-02',
                kn.mrn_col: '01',
                kn.vital_status_col: 'alive'
            }
        ]
        me.create_trial_match_records()

        trial_matches = self._findall(query={},
                                      proj={k: 1 for k in trial_matches_schema.keys()},
                                      table=s.trial_match_collection_name)
        self._print(trial_matches)
        for tm in trial_matches:
            assert tm[kn.tm_trial_protocol_no_col] == '00-000'
            assert tm[kn.sample_id_col] == 'DEV-01' or tm[kn.sample_id_col] == 'DEV-02'
            assert tm[kn.mrn_col] == '01'
            assert tm[kn.tm_sort_order_col] == 0
            assert tm[kn.vital_status_col] == 'alive'
            assert tm[kn.tm_trial_accrual_status_col] == 'open'
            assert tm[kn.mr_trial_level_col] == 'step'
            assert tm[kn.mr_trial_step_code_col] == 'STEP CODE'
            assert tm[kn.mr_trial_arm_code_col] is None
            assert tm[kn.mr_trial_dose_code_col] is None
            assert kn.mutation_list_col in tm
            assert tm[kn.oncotree_primary_diagnosis_name_col] == 'Lung Adenocarcinoma'
            assert tm[kn.birth_date_col] == dt.datetime(year=1900, day=1, month=1)
            assert tm[kn.mr_diagnosis_level_col] == '_solid_' or tm[kn.mr_diagnosis_level_col] == 'specific'
