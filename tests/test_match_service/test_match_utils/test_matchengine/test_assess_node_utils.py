from src.utilities import settings as s
from src.data_store import key_names as kn
from tests.test_match_service import TestQueryUtilitiesShared
from src.services.match_service.match_utils.matchengine.assess_node_utils import AssessNodeUtils


class TestAssessNodeUtils(TestQueryUtilitiesShared):

    def setUp(self):
        super(TestAssessNodeUtils, self).setUp()
        self.a = AssessNodeUtils()
        self.query = {'$and': []}
        self.proj_info = []

        self.standard_genomic_mut_proj = {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.cnv_call_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.sv_comment_col): 1,
        }
        self.standard_genomic_cnv_proj = {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            '%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.protein_change_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.variant_class_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.cnv_call_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.ref_residue_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.transcript_exon_col): 1,
            '%s.%s' % (kn.cnv_list_col, kn.sv_comment_col): 1,
        }
        self.standard_genomic_sv_proj = {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            '%s.%s' % (kn.sv_list_col, kn.hugo_symbol_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.protein_change_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.variant_class_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.cnv_call_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.ref_residue_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.transcript_exon_col): 1,
            '%s.%s' % (kn.sv_list_col, kn.sv_comment_col): 1,
        }
        self.standard_genomic_wt_proj = {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            '%s.%s' % (kn.wt_genes_col, kn.hugo_symbol_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.protein_change_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.variant_class_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.cnv_call_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.ref_residue_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.transcript_exon_col): 1,
            '%s.%s' % (kn.wt_genes_col, kn.sv_comment_col): 1,
        }

    def tearDown(self):
        pass

    def test_parse_diagnosis(self):

        # inclusion
        node = {'value': {s.mt_diagnosis: 'Lung Adenocarcinoma'}}
        i = self.a._parse_diagnosis(node=node, query=self.query, proj_info=self.proj_info)
        assert i is True
        assert self.query['$and'][0] == {kn.oncotree_primary_diagnosis_name_col: {'$in': ['Lung Adenocarcinoma']}}
        assert self.proj_info[0] == {s.mt_diagnosis: 'Lung Adenocarcinoma'}
        assert node[kn.mr_diagnosis_level_col] == 'specific'

        # exclusion
        node = {'value': {s.mt_diagnosis: '!Lung Adenocarcinoma'}}
        i = self.a._parse_diagnosis(node=node, query=self.query, proj_info=self.proj_info)
        assert i is False
        assert self.query['$and'][1] == {kn.oncotree_primary_diagnosis_name_col: {'$nin': ['Lung Adenocarcinoma']}}
        assert self.proj_info[1] == {s.mt_diagnosis: '!Lung Adenocarcinoma'}
        assert node[kn.mr_diagnosis_level_col] == 'specific'

        # _solid_ and _liquid_ expansion capture
        node = {'value': {s.mt_diagnosis: s.oncotree_all_solid_text}}
        i = self.a._parse_diagnosis(node=node, query=self.query, proj_info=self.proj_info)
        assert i is True
        assert node[kn.mr_diagnosis_level_col] == '_solid_'

        node = {'value': {s.mt_diagnosis: s.oncotree_all_liquid_text}}
        i = self.a._parse_diagnosis(node=node, query=self.query, proj_info=self.proj_info)
        assert i is True
        assert node[kn.mr_diagnosis_level_col] == '_liquid_'

    def test_parse_age(self):

        node = {'value': {s.mt_age: '>=18'}}
        self.a._parse_age(node=node, query=self.query, proj_info=self.proj_info)
        assert self.query['$and'][0].keys() == [kn.birth_date_col]
        assert self.query['$and'][0][kn.birth_date_col].keys() == ['$lte']
        assert self.proj_info[0] == {s.mt_age: '>=18'}

    def test_parse_gender(self):

        node = {'value': {s.mt_gender: 'Female'}}
        self.a._parse_gender(node=node, query=self.query, proj_info=self.proj_info)
        assert self.query == {'$and': [{kn.gender_col: 'Female'}]}
        assert self.proj_info == [{s.mt_gender: 'Female'}]

    def test_parse_variant_category(self):

        # inclusive mutation
        node = {'value': {s.mt_variant_category: s.mt_mut_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_mutation_val
        assert i is True

        # exclusive mutation
        node = {'value': {s.mt_variant_category: '!%s' % s.mt_mut_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_mutation_val
        assert i is False

        # inclusive cnv
        node = {'value': {s.mt_variant_category: s.mt_cnv_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_cnv_val
        assert i is True

        # exclusive cnv
        node = {'value': {s.mt_variant_category: '!%s' % s.mt_cnv_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_cnv_val
        assert i is False

        # inclusive sv
        node = {'value': {s.mt_variant_category: s.mt_sv_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_sv_val
        assert i is True

        # exclusive sv
        node = {'value': {s.mt_variant_category: '!%s' % s.mt_sv_val}}
        vc, i = self.a._parse_variant_category(node=node)
        assert vc == s.variant_category_sv_val
        assert i is False

    def test_parse_gene_level(self):

        # inclusive mutation
        node = {'value': {
            s.mt_variant_category: s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_mut_proj
        assert node['variant_level'] == 'gene'
        assert node['unwind'] == '$%s' % kn.mutation_list_col
        assert node['match_reason'] == {'%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF'}
        assert node['include'] is True

        # exclusive mutation
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_mutation_val,
            kn.hugo_symbol_col: 'BRAF'
        }
        assert node['variant_level'] == 'gene'
        assert node['include'] is False

        # inclusive cnv
        node = {'value': {
            s.mt_variant_category: s.mt_cnv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_cnv_proj
        assert node['variant_level'] == 'gene'
        assert node['unwind'] == '$%s' % kn.cnv_list_col
        assert node['match_reason'] == {'%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF'}
        assert node['include'] is True

        # exclusive cnv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_cnv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_cnv_val,
            kn.hugo_symbol_col: 'BRAF'
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'gene'
        assert node['include'] is False

        # inclusive sv
        node = {'value': {
            s.mt_variant_category: s.mt_sv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_sv_proj
        assert node['variant_level'] == 'gene'
        assert node['unwind'] == '$%s' % kn.sv_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.sv_list_col, kn.sv_comment_col): node['query'][kn.sv_list_col]['$elemMatch'][kn.sv_comment_col]
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_sv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_sv_val,
            kn.hugo_symbol_col: 'BRAF'
        }
        assert node['variant_level'] == 'gene'
        assert node['include'] is False

        # todo enable
        # # inclusive any variant
        # node = {'value': {
        #     s.mt_variant_category: s.mt_any_vc_val,
        #     s.mt_hugo_symbol: 'BRAF'}
        # }
        # self.a._parse_gene_level(node=node)
        # assert 'query' in node
        # assert 'genomic_inclusion_reasons' in node
        # print node['genomic_inclusion_reasons']
        # assert node['genomic_inclusion_reasons'] == {
        #     '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.cnv_call_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 1,
        #     '%s.%s' % (kn.mutation_list_col, kn.sv_comment_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.protein_change_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.variant_class_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.cnv_call_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.ref_residue_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.sv_comment_col): 1,
        #     '%s.%s' % (kn.cnv_list_col, kn.transcript_exon_col): 1,
        # }
        # assert node['variant_level'] == 'gene'
        # assert node['unwind'] == '$%s' % kn.mutation_list_col  # todo bad
        # assert node['match_reason'] == {
        #     '$or': [
        #         {'%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF'},
        #         {'%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF'}
        #     ]
        # }
        #
        # # exclusive any variant
        # node = {'value': {
        #     s.mt_variant_category: '!%s' % s.mt_any_vc_val,
        #     s.mt_hugo_symbol: 'BRAF'}
        # }
        # self.a._parse_gene_level(node=node)
        # assert 'query' in node
        # assert kn.genomic_exclusion_reasons_col in node
        # assert node[kn.genomic_exclusion_reasons_col] == {
        #     '$and': [
        #         {kn.variant_category_col: s.variant_category_mutation_val, kn.hugo_symbol_col: 'BRAF'},
        #         {kn.variant_category_col: s.variant_category_cnv_val, kn.hugo_symbol_col: 'BRAF'}
        #     ]
        # }
        # assert node['variant_level'] == 'gene'

    def test_parse_sv(self):

        # inclusive sv
        node = {'value': {
            s.mt_variant_category: s.mt_sv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_sv_proj
        assert node['variant_level'] == 'gene'
        assert node['unwind'] == '$%s' % kn.sv_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.sv_list_col, kn.sv_comment_col): node['query'][kn.sv_list_col]['$elemMatch'][
                kn.sv_comment_col]
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_sv_val,
            s.mt_hugo_symbol: 'BRAF'}
        }
        self.a._parse_gene_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_sv_val,
            kn.hugo_symbol_col: 'BRAF'
        }
        assert node['variant_level'] == 'gene'
        assert node['include'] is False

    def test_parse_variant_level(self):

        # inclusive
        node = {'value': {
            s.mt_variant_category: s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_protein_change: 'p.V600E'}
        }
        self.a._parse_variant_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_mut_proj
        assert node['variant_level'] == 'variant'
        assert node['unwind'] == '$%s' % kn.mutation_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 'p.V600E'
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_protein_change: 'p.V600E'}
        }
        self.a._parse_variant_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_mutation_val,
            kn.hugo_symbol_col: 'BRAF',
            kn.protein_change_col: 'p.V600E'
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'variant'
        assert node['include'] is False

    def test_parse_wildcard_level(self):

        # inclusive
        node = {'value': {
            s.mt_variant_category: s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_wc_protein_change: 'p.V600'}
        }
        self.a._parse_wildcard_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_mut_proj
        assert node['variant_level'] == 'wildcard'
        assert node['unwind'] == '$%s' % kn.mutation_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 'p.V600'
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_wc_protein_change: 'p.V600'}
        }
        self.a._parse_wildcard_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_mutation_val,
            kn.hugo_symbol_col: 'BRAF',
            kn.ref_residue_col: 'p.V600'
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'wildcard'
        assert node['include'] is False

    def test_parse_exon_level(self):

        # inclusive
        node = {'value': {
            s.mt_variant_category: s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_exon: 20}
        }
        self.a._parse_exon_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_mut_proj
        assert node['variant_level'] == 'exon'
        assert node['unwind'] == '$%s' % kn.mutation_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 20
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_exon: 20}
        }
        self.a._parse_exon_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_mutation_val,
            kn.hugo_symbol_col: 'BRAF',
            kn.transcript_exon_col: 20
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'exon'
        assert node['include'] is False

    def test_parse_variant_class_level(self):

        # inclusive
        node = {'value': {
            s.mt_variant_category: s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_variant_class: 'Nonsense_Mutation'}
        }
        self.a._parse_variant_class_level(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_mut_proj
        assert node['variant_level'] == 'variant_class'
        assert node['unwind'] == '$%s' % kn.mutation_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 'Nonsense_Mutation'
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_mut_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_variant_class: 'Nonsense_Mutation'}
        }
        self.a._parse_variant_class_level(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_mutation_val,
            kn.hugo_symbol_col: 'BRAF',
            kn.variant_class_col: 'Nonsense_Mutation'
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'variant_class'
        assert node['include'] is False

    def test_parse_cnv_call(self):

        # inclusive
        node = {'value': {
            s.mt_variant_category: s.mt_cnv_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_cnv_call: s.mt_hetero_del_val}
        }
        self.a._parse_cnv_call(node=node)
        assert 'query' in node
        assert 'genomic_inclusion_reasons' in node
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_cnv_proj
        assert node['variant_level'] == 'variant'
        assert node['unwind'] == '$%s' % kn.cnv_list_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.cnv_list_col, kn.cnv_call_col): s.cnv_call_hetero_del
        }
        assert node['include'] is True

        # exclusive sv
        node = {'value': {
            s.mt_variant_category: '!%s' % s.mt_cnv_val,
            s.mt_hugo_symbol: 'BRAF',
            s.mt_cnv_call: s.mt_hetero_del_val}
        }
        self.a._parse_cnv_call(node=node)
        assert 'query' in node
        assert kn.genomic_exclusion_reasons_col in node
        assert node[kn.genomic_exclusion_reasons_col] == {
            kn.variant_category_col: s.variant_category_cnv_val,
            kn.hugo_symbol_col: 'BRAF',
            kn.cnv_call_col: s.cnv_call_hetero_del
        }, node[kn.genomic_exclusion_reasons_col]
        assert node['variant_level'] == 'variant'
        assert node['include'] is False

    def test_parse_signature(self):

        # MMR Status
        node = {'value': {s.mt_mmr_status: s.mt_mmr_deficient_val}}
        self.a._parse_signature(node=node, criteria=[s.mt_mmr_status])
        assert node['query'] == {kn.mmr_status_col: s.mmr_status_deficient_val}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.mmr_status_col: 1
        }
        assert node['match_reason'][kn.mmr_status_col] == s.mmr_status_deficient_val

        # MS Status
        node = {'value': {s.mt_ms_status: s.mt_msi_high_val}}
        self.a._parse_signature(node=node, criteria=[s.mt_ms_status])
        assert node['query'] == {kn.ms_status_col: s.ms_status_msih_val}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.ms_status_col: 1
        }
        assert node['match_reason'][kn.ms_status_col] == s.mt_msi_high_val

        # Tobacco Status
        node = {'value': {s.mt_tobacco_status: 'Yes'}}
        self.a._parse_signature(node=node, criteria=[s.mt_tobacco_status])
        assert node['query'] == {kn.tobacco_status_col: 'Yes'}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.tobacco_status_col: 1
        }
        assert node['match_reason'][kn.tobacco_status_col] == 'Yes'

        # TMZ Status
        node = {'value': {s.mt_tmz_status: 'Yes'}}
        self.a._parse_signature(node=node, criteria=[s.mt_tmz_status])
        assert node['query'] == {kn.tmz_status_col: 'Yes'}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.tmz_status_col: 1
        }
        assert node['match_reason'][kn.tmz_status_col] == 'Yes'

        # polE Status
        node = {'value': {s.mt_pole_status: 'Yes'}}
        self.a._parse_signature(node=node, criteria=[s.mt_pole_status])
        assert node['query'] == {kn.pole_status_col: 'Yes'}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.pole_status_col: 1
        }
        assert node['match_reason'][kn.pole_status_col] == 'Yes'

        # APOBEC Status
        node = {'value': {s.mt_apobec_status: 'Yes'}}
        self.a._parse_signature(node=node, criteria=[s.mt_apobec_status])
        assert node['query'] == {kn.apobec_status_col: 'Yes'}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.apobec_status_col: 1
        }
        assert node['match_reason'][kn.apobec_status_col] == 'Yes'

        # UVA Status
        node = {'value': {s.mt_uva_status: 'Yes'}}
        self.a._parse_signature(node=node, criteria=[s.mt_uva_status])
        assert node['query'] == {kn.uva_status_col: 'Yes'}
        assert node['genomic_inclusion_reasons'] == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.uva_status_col: 1
        }
        assert node['match_reason'][kn.uva_status_col] == 'Yes'

    def test_parse_wildtype(self):

        node = {'value': {s.mt_hugo_symbol: 'BRAF', s.mt_wildtype: True}}
        self.a._parse_wildtype(node=node)
        assert node['query'] == {kn.wt_genes_col: {'$elemMatch': {kn.hugo_symbol_col: 'BRAF'}}}
        assert node['genomic_inclusion_reasons'] == self.standard_genomic_wt_proj
        assert node['unwind'] == '$%s' % kn.wt_genes_col
        assert node['match_reason'] == {
            '%s.%s' % (kn.wt_genes_col, kn.hugo_symbol_col): 'BRAF'
        }
        assert node['include'] is True
