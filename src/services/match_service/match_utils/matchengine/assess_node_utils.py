from src.utilities import settings as s
from src.data_store import key_names as kn
from src.services.match_service.match_utils.matchengine import matchengine_utils as me_utils
from src.services.match_service.query_utils.clinical_queries import ClinicalQueries
from src.services.match_service.query_utils.genomic_queries import GenomicQueries
from src.services.match_service.query_utils.proj_utils import ProjUtils


class AssessNodeUtils(ClinicalQueries, GenomicQueries, ProjUtils):

    def __init__(self):
        ClinicalQueries.__init__(self)
        GenomicQueries.__init__(self)
        ProjUtils.__init__(self)

        self.proj_dict = {True: 'inclusion_reasons', False: 'exclusion_reasons'}

    def assess_clinical_node(self, node):
        """
        Assess the given node and construct the appropriate MongoDB query.

        :param node: {digraph node}
        :return: {digraph node}
        """
        query = {'$and': []}
        proj_info = []
        criteria = sorted(node['value'].keys())

        if s.mt_diagnosis not in node['value']:
            raise ValueError('%s column must be included for all clinical nodes' % s.mt_diagnosis)

        # diagnosis query
        include = True
        if s.mt_diagnosis in criteria:
            include = self._parse_diagnosis(node=node, query=query, proj_info=proj_info)

        # age query
        if s.mt_age in criteria:
            self._parse_age(node=node, query=query, proj_info=proj_info)

        # gender query
        if s.mt_gender in criteria:
            self._parse_gender(node=node, query=query, proj_info=proj_info)

        # clinical projection
        proj_keys = [i.keys()[0] for i in proj_info]
        proj_vals = [i.values()[0] for i in proj_info]
        node['clinical_%s' % self.proj_dict[include]] = self.create_clinical_proj(include=include,
                                                                                  keys=proj_keys,
                                                                                  vals=proj_vals)

        node['query'] = query
        return node

    def assess_genomic_node(self, node):
        """
        Assess the given node and construct the appropriate MongoDB query.

        :param node: {digraph node}
        :return: {digraph node}
        """
        ignore_these_fields = ['display_name', 'fusion_partner_hugo_symbol']
        criteria = sorted([i for i in node['value'].keys() if i not in ignore_these_fields])

        # gene level query (mutations, cnvs, and svs)
        if criteria == [s.mt_hugo_symbol, s.mt_variant_category]:
            return self._parse_gene_level(node=node)

        # variant-level mutation criteria
        elif s.mt_protein_change in criteria:
            return self._parse_variant_level(node=node)

        # wildcard-level mutation criteria
        elif s.mt_wc_protein_change in criteria:
            return self._parse_wildcard_level(node=node)

        # exon-level mutation criteria
        elif s.mt_exon in criteria:
            return self._parse_exon_level(node=node)

        # variant-class-level mutation criteria
        elif s.mt_variant_class in criteria:
            return self._parse_variant_class_level(node=node)

        # cnv criteria
        elif s.mt_cnv_call in criteria:
            return self._parse_cnv_call(node=node)

        # mutational signature criteria
        elif any([criterion in s.mt_signature_cols for criterion in criteria]):
            return self._parse_signature(node=node, criteria=criteria)

        # wildtype criteria
        elif s.mt_wildtype in node['value'] and node['value'][s.mt_wildtype] is True:
            return self._parse_wildtype(node=node)

        # low-coverage criteria
        # todo build low coverage criteria logic
        else:
            raise ValueError('This node does not match an expected format.')

    def _parse_diagnosis(self, node, query, proj_info):
        """
        Parse cancer type text from node and add to the query and projection info.

        :param node {digraph node}
        :param query {dict}
        :param proj_info {dict}
        :return: {bool}
        """
        cancer_type = node['value'][s.mt_diagnosis]
        include = me_utils.assess_inclusion(cancer_type)
        sanitized_cancer_type = me_utils.sanitize_exclusion_vals(cancer_type)
        subquery = self.create_oncotree_diagnosis_query(cancer_type=sanitized_cancer_type,
                                                        include=include)
        node[kn.mr_diagnosis_level_col] = sanitized_cancer_type.lower() \
            if cancer_type in [s.oncotree_all_solid_text, s.oncotree_all_liquid_text] \
            else 'specific'

        query['$and'].append(subquery)
        proj_info.append({s.mt_diagnosis: cancer_type})
        return include

    def _parse_age(self, node, query, proj_info):
        """
        Parse age text from node and add to the query and projection info.

        :param node: {digraph node}
        :param query: {dict}
        :param proj_info: {dict}
        :return: {null}
        """
        age = node['value'][s.mt_age]
        subquery = self.create_age_query(age=age)
        query['$and'].append(subquery)
        proj_info.append({s.mt_age: age})

    def _parse_gender(self, node, query, proj_info):
        """
        Parse gender text from node and add to the query and projection info.

        :param node: {digraph node}
        :param query: {dict}
        :param proj_info: {dict}
        :return: {null}
        """
        gender = node['value'][s.mt_gender]
        subquery = self.create_gender_query(gender=gender)
        query['$and'].append(subquery)
        proj_info.append({s.mt_gender: gender})

    @staticmethod
    def _parse_variant_category(node):
        """
        Parse variant category and include info from node

        :param node: {digraph node}
        :return: {tuple} (variant_category {str}, include {bool})
        """
        include = True
        variant_category = None
        if s.mt_variant_category in node['value']:
            variant_category = me_utils.sanitize_exclusion_vals(node['value'][s.mt_variant_category])
            variant_category = me_utils.normalize_variant_category_val(variant_category)
            include = me_utils.assess_inclusion(node_value=node['value'][s.mt_variant_category])

        return variant_category, include

    def _parse_gene_level(self, node):
        """
        Parse gene-level mutation, cnv, and sv criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'gene'
        gene_name = node['value'][s.mt_hugo_symbol]
        variant_category, include = self._parse_variant_category(node=node)
        if variant_category == s.variant_category_sv_val:
            return me_utils.add_unwind(node=self._parse_sv(node=node), include=include)

        proj = 'genomic_%s' % self.proj_dict[include]

        # Any Variation
        if variant_category == s.variant_category_any_val:

            # query
            node['query'] = self.create_any_variant_query(gene_name=gene_name, include=include)

            # projection
            node[proj], node['match_reason'] = self.create_any_variant_proj(gene_name=gene_name,
                                                                            include=include,
                                                                            query=node['query'])
            node['unwind'] = '$%s' % kn.mutation_list_col

        else:

            # query
            node['query'] = self.create_gene_level_query(gene_name=gene_name,
                                                         variant_category=variant_category,
                                                         include=include)

            # projection
            proj = 'genomic_%s' % self.proj_dict[include]
            proj_info = {
                self.variant_category_dict[variant_category]: variant_category,
                self.hugo_symbol_key: gene_name
            }
            node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                        query=node['query'],
                                                                        keys=proj_info.keys(),
                                                                        vals=proj_info.values())

        return me_utils.add_unwind(node=node, include=include)

    def _parse_sv(self, node):
        """
        Parse structural variant criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        gene_name = node['value'][s.mt_hugo_symbol]
        variant_category, include = self._parse_variant_category(node=node)

        # query
        node['query'] = self.create_sv_query(gene_name=gene_name, include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name
        }
        node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                    query=node['query'],
                                                                    keys=proj_info.keys(),
                                                                    vals=proj_info.values())
        return node

    def _parse_variant_level(self, node):
        """
        Parse variant-level mutation criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'variant'
        gene_name = node['value'][s.mt_hugo_symbol]
        protein_change = node['value'][s.mt_protein_change]
        variant_category, include = self._parse_variant_category(node=node)

        # query
        node['query'] = self.create_mutation_query(gene_name=gene_name,
                                                   protein_change=protein_change,
                                                   include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name,
            self.protein_change_key: protein_change
        }
        node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                    query=node['query'],
                                                                    keys=proj_info.keys(),
                                                                    vals=proj_info.values())
        return me_utils.add_unwind(node=node, include=include)

    def _parse_wildcard_level(self, node):
        """
        Parse wildcard-level mutation criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'wildcard'
        gene_name = node['value'][s.mt_hugo_symbol]
        protein_change = node['value'][s.mt_wc_protein_change]
        variant_category, include = self._parse_variant_category(node=node)

        # query
        node['query'] = self.create_wildcard_query(gene_name=gene_name,
                                                   protein_change=protein_change,
                                                   include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name,
            self.ref_residue_key: protein_change
        }
        node[proj],  node['match_reason'] = self.create_genomic_proj(include=include,
                                                                     query=node['query'],
                                                                     keys=proj_info.keys(),
                                                                     vals=proj_info.values())

        return me_utils.add_unwind(node=node, include=include)

    def _parse_exon_level(self, node):
        """
        Parse exon-level mutation criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'exon'
        exon = node['value'][s.mt_exon]
        gene_name = node['value'][s.mt_hugo_symbol]
        variant_category, include = self._parse_variant_category(node=node)
        variant_class = node['value'][s.mt_variant_class] if s.mt_variant_class in node['value'] else None

        # query
        node['query'] = self.create_exon_query(gene_name=gene_name,
                                               exon=exon,
                                               variant_class=variant_class,
                                               include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name,
            self.transcript_exon_key: exon,
            self.variant_class_key: variant_class
        }
        node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                    query=node['query'],
                                                                    keys=proj_info.keys(),
                                                                    vals=proj_info.values())

        return me_utils.add_unwind(node=node, include=include)

    def _parse_variant_class_level(self, node):
        """
        Parse variant-class-level mutation criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'variant_class'
        variant_class = node['value'][s.mt_variant_class]
        gene_name = node['value'][s.mt_hugo_symbol]
        variant_category, include = self._parse_variant_category(node=node)

        # query
        node['query'] = self.create_variant_class_query(gene_name=gene_name,
                                                        variant_class=variant_class,
                                                        include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name,
            self.variant_class_key: variant_class
        }
        node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                    query=node['query'],
                                                                    keys=proj_info.keys(),
                                                                    vals=proj_info.values())

        return me_utils.add_unwind(node=node, include=include)

    def _parse_cnv_call(self, node):
        """
        Parse cnv call criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        node['variant_level'] = 'variant'
        gene_name = node['value'][s.mt_hugo_symbol]
        variant_category, include = self._parse_variant_category(node=node)
        cnv_call = me_utils.normalize_cnv_call_val(me_utils.sanitize_exclusion_vals(node['value'][s.mt_cnv_call]))

        # query
        node['query'] = self.create_cnv_query(gene_name=gene_name,
                                              cnv_call=cnv_call,
                                              include=include)

        # projection
        proj = 'genomic_%s' % self.proj_dict[include]
        proj_info = {
            self.variant_category_dict[variant_category]: variant_category,
            self.hugo_symbol_key: gene_name,
            self.cnv_call_key: cnv_call
        }
        node[proj], node['match_reason'] = self.create_genomic_proj(include=include,
                                                                    query=node['query'],
                                                                    keys=proj_info.keys(),
                                                                    vals=proj_info.values())

        return me_utils.add_unwind(node=node, include=include)

    def _parse_signature(self, node, criteria):
        """
        Parse signature criteria from node

        :param node: {digraph node}
        :param criteria: {list of str}
        :return: {digraph node}
        """
        for sig in s.mt_signature_cols:
            if sig in criteria:
                sigtype, sigval = me_utils.normalize_signature_vals(signature_type=sig, signature_val=node['value'][sig])
                node['query'] = self.create_mutational_signature_query(signature_type=sigtype, signature_val=sigval)
                node['genomic_inclusion_reasons'], node['match_reason'] = self.create_genomic_proj(include=True,
                                                                                                   query=node['query'])
                return node

    def _parse_wildtype(self, node):
        """
        Parse wildtype criteria from node

        :param node: {digraph node}
        :return: {digraph node}
        """
        # parse node
        gene_name = node['value'][s.mt_hugo_symbol]

        # query
        node['query'] =self.create_gene_level_query(gene_name=gene_name,
                                                    variant_category=s.variant_category_wt_val,
                                                    include=True)
        # projection
        node['genomic_inclusion_reasons'], node['match_reason'] = self.create_genomic_proj(include=True,
                                                                                           query=node['query'])
        return me_utils.add_unwind(node=node, include=True)
