"""Copyright 2016 Dana-Farber Cancer Institute"""

from tests.test_match_service import TestQueryUtilitiesShared, get_demo_trial_matches
from src.services.match_service.match_utils.sort import *


class TestSort(TestQueryUtilitiesShared):

    def setUp(self):
        super(TestSort, self).setUp()

        self.s = Sort(trial_matches=[])

    def tearDown(self):
        pass

    def test_has_vc(self):

        assert has_vc(None) is False
        assert has_vc([]) is False
        assert has_vc([[1]]) is True

    def test_prep_trial_matches(self):

        data = [
            {
                kn.sample_id_col: 'DEV-01',
                kn.vital_status_col: 'alive',
                kn.tm_trial_accrual_status_col: 'open',
                kn.mutation_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
                kn.cnv_list_col: [],
                kn.sv_list_col: [],
                kn.wt_genes_col: None
            },
            {
                kn.sample_id_col: 'DEV-02',
                kn.vital_status_col: 'alive',
                kn.tm_trial_accrual_status_col: 'open',
                kn.mutation_list_col: [],
                kn.cnv_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
                kn.sv_list_col: [],
                kn.wt_genes_col: None
            },
            {
                kn.sample_id_col: 'DEV-03',
                kn.vital_status_col: 'alive',
                kn.tm_trial_accrual_status_col: 'open',
                kn.mutation_list_col: [],
                kn.cnv_list_col: [],
                kn.sv_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
                kn.wt_genes_col: None
            },
            {
                kn.sample_id_col: 'DEV-03',
                kn.mutation_list_col: [],
                kn.cnv_list_col: [],
                kn.sv_list_col: None,
                kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}]
            }
        ]
        s = Sort(trial_matches=data)
        s._prep_trial_matches()
        assert s.trial_matches_df['has_mut'].tolist() == [True, False, False, False]
        assert s.trial_matches_df['has_cnv'].tolist() == [False, True, False, False]
        assert s.trial_matches_df['has_sv'].tolist() == [False, False, True, False]
        assert s.trial_matches_df['has_wt'].tolist() == [False, False, False, True]
        assert sorted(s.all_sample_ids) == ['DEV-01', 'DEV-02', 'DEV-03']
        assert s.f_alive is not None
        assert s.f_open is not None
        assert s.f_no_svs is not None

    def test_add_sort_value(self):

        li = []
        li = add_sort_value(8, 0, li)
        assert li == [8]

        li = add_sort_value(0, 0, li)
        assert li == [0]

        li = add_sort_value(1, 1, li)
        assert li == [0, 1]

        li = [0, 1, 0, 0]
        li = add_sort_value(1, 0, li)
        assert li == [0, 1, 0, 0]

    def test_extract_tier(self):

        m1 = {}
        m2 = {kn.mutation_list_col: []}
        m3 = {kn.mutation_list_col: [{}]}
        m4 = {kn.mutation_list_col: [{kn.tier_col: 1}]}
        m5 = {kn.mutation_list_col: [{kn.tier_col: 1}, {}]}
        m6 = {kn.mutation_list_col: [{kn.tier_col: 1}, {kn.tier_col: None}]}
        m7 = {kn.mutation_list_col: [{kn.tier_col: 1}, {kn.tier_col: 2}]}
        assert extract_tier(m1) is None
        assert extract_tier(m2) is None
        assert extract_tier(m3) is None
        assert extract_tier(m4) == 1
        assert extract_tier(m5) == 1
        assert extract_tier(m6) == 1
        assert extract_tier(m7) == 1

    def test_extract_match_level(self):

        m1 = {}
        m2 = {kn.mutation_list_col: []}
        m3 = {kn.mutation_list_col: [{}]}
        m4 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'variant'}]}
        m5 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'variant'}, {}]}
        m6 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'variant'}, {kn.mr_reason_level_col: 'wildcard'}]}
        m7 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'variant'}, {kn.mr_reason_level_col: 'wildcard'},
                                     {kn.mr_reason_level_col: 'exon'}, {kn.mr_reason_level_col: 'gene'}]}
        m8 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'wildcard'}, {kn.mr_reason_level_col: 'exon'},
                                     {kn.mr_reason_level_col: 'gene'}]}
        m9 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'exon'}, {kn.mr_reason_level_col: 'gene'}]}
        m10 = {kn.mutation_list_col: [{kn.mr_reason_level_col: 'gene'}]}
        assert extract_match_level(m1) is None
        assert extract_match_level(m2) is None
        assert extract_match_level(m3) is None
        assert extract_match_level(m4) == 'variant'
        assert extract_match_level(m5) == 'variant'
        assert extract_match_level(m6) == 'variant'
        assert extract_match_level(m7) == 'variant'
        assert extract_match_level(m8) == 'wildcard'
        assert extract_match_level(m9) == 'exon'
        assert extract_match_level(m10) == 'gene'

    def test_sort_by_tier(self):

        m1 = ('DEV-01', '00-001')

        # MMR Deficient
        sort_order = {m1: []}
        m1_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mmr_status_col: s.mmr_status_deficient_val
        }
        sort_order = sort_by_tier(m1_data, sort_order)
        assert sort_order[m1][0] == 0

        # Tier 1 BRAF p.V600E
        sort_order = {m1: []}
        m2_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.cnv_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.V600E', kn.tier_col: 1},
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.G469A', kn.tier_col: 2}
            ]
        }
        sort_order = sort_by_tier(m2_data, sort_order)
        assert sort_order[m1][0] == 1

        # Tier 2 BRAF p.G469A
        sort_order = {m1: []}
        m3_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.cnv_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.G469A', kn.tier_col: 2}
            ]
        }
        sort_order = sort_by_tier(m3_data, sort_order)
        assert sort_order[m1][0] == 2

        # Tier 3 BRAF p.N581S
        m4_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.N581S', kn.tier_col: 3}
            ]
        }
        sort_order = {m1: []}
        sort_order = sort_by_tier(m4_data, sort_order)
        assert sort_order[m1][0] == 4

        # Tier 4 BRAF p.E695Q
        sort_order = {m1: []}
        m5_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.E695Q', kn.tier_col: 4}
            ]
        }
        sort_order = sort_by_tier(m5_data, sort_order)
        assert sort_order[m1][0] == 5

        # Weird unknown edge case
        sort_order = {m1: []}
        m6_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
        }
        sort_order = sort_by_tier(m6_data, sort_order)
        assert sort_order[m1][0] == 7

        # BRAF CNV
        sort_order = {m1: []}
        m7_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.cnv_list_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}],
            kn.mutation_list_col: [
                {kn.hugo_symbol_col: 'BRAF', kn.protein_change_col: 'p.E695Q', kn.tier_col: 4}
            ]
        }
        sort_order = sort_by_tier(m7_data, sort_order)
        assert sort_order[m1][0] == 3

        # Wild-type BRAF
        sort_order = {m1: []}
        m8_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.wt_genes_col: [{kn.hugo_symbol_col: 'BRAF'}]
        }
        sort_order = sort_by_tier(m8_data, sort_order)
        assert sort_order[m1][0] == 6

    def test_sort_by_match_type(self):

        m1 = ('DEV-01', '00-001')

        # variant-level
        sort_order = {m1: [0]}
        m1_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mutation_list_col: [{kn.mr_reason_level_col: 'variant'}]
        }
        sort_order = sort_by_match_type(m1_data, sort_order)
        assert sort_order[m1][1] == 0

        # wildcard-level
        sort_order = {m1: [0]}
        m2_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mutation_list_col: [{kn.mr_reason_level_col: 'wildcard'}]
        }
        sort_order = sort_by_match_type(m2_data, sort_order)
        assert sort_order[m1][1] == 1

        # exon-level
        sort_order = {m1: [0]}
        m3_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mutation_list_col: [{kn.mr_reason_level_col: 'exon'}]
        }
        sort_order = sort_by_match_type(m3_data, sort_order)
        assert sort_order[m1][1] == 2

        # gene-level
        sort_order = {m1: [0]}
        m4_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mutation_list_col: [{kn.mr_reason_level_col: 'gene'}]
        }
        sort_order = sort_by_match_type(m4_data, sort_order)
        assert sort_order[m1][1] == 3

        # no mutation
        sort_order = {m1: [0]}
        m5_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
        }
        sort_order = sort_by_match_type(m5_data, sort_order)
        assert sort_order[m1][1] == 4

    def test_sort_by_cancer_type(self):

        m1 = ('DEV-01', '00-001')

        # specifically Lung Adenocarcinoma cancer type
        sort_order = {m1: [0, 0]}
        m1_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.mr_diagnosis_level_col: 'specific'
        }
        sort_order = sort_by_cancer_type(m1_data, sort_order)
        assert sort_order[m1][2] == 0

        # all sold tumor expansion match
        sort_order = {m1: [0, 0]}
        m2_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.mr_diagnosis_level_col: '_solid_'
        }
        sort_order = sort_by_cancer_type(m2_data, sort_order)
        assert sort_order[m1][2] == 1

        # all liquid tumor expansion match
        sort_order = {m1: [0, 0]}
        m3_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
            kn.mr_diagnosis_level_col: '_liquid_'
        }
        sort_order = sort_by_cancer_type(m3_data, sort_order)
        assert sort_order[m1][2] == 1

        # cancer type level not specified
        sort_order = {m1: [0, 0]}
        m4_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.oncotree_primary_diagnosis_name_col: 'Lung Adenocarcinoma',
        }
        sort_order = sort_by_cancer_type(m4_data, sort_order)
        assert sort_order[m1][2] == 2

    def test_sort_by_coordinating_center(self):

        m1 = ('DEV-01', '00-001')

        # DFCI
        sort_order = {m1: [0, 0]}
        m1_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mr_coordinating_center_col: s.coordinating_center_dfci_val
        }
        sort_order = sort_by_coordinating_center(m1_data, sort_order)
        assert sort_order[m1][2] == 0

        # MGH
        sort_order = {m1: [0, 0]}
        m2_data = {
            kn.mrn_col: 'MRN-01',
            kn.sample_id_col: m1[0],
            kn.tm_trial_protocol_no_col: m1[1],
            kn.mr_coordinating_center_col: 'Massachussetts General Hospital'
        }
        sort_order = sort_by_coordinating_center(m2_data, sort_order)
        assert sort_order[m1][2] == 1

    def test_sort_by_reverse_protocol_no(self):

        sort_order = {
            ('DEV-01', '11-111'): [0, 0, 0, 0],
            ('DEV-01', '09-999'): [0, 0, 0, 1],
            ('DEV-01', '15-000'): [7, 0, 0, 0],
            ('DEV-01', '15-111'): [7, 0, 0, 0],
            ('DEV-01', '22-222'): [7, 1, 0, 0],
        }
        matches = [
            {kn.tm_trial_protocol_no_col: '11-111', kn.sample_id_col: 'DEV-01'},
            {kn.tm_trial_protocol_no_col: '09-999', kn.sample_id_col: 'DEV-01'},
            {kn.tm_trial_protocol_no_col: '15-000', kn.sample_id_col: 'DEV-01'},
            {kn.tm_trial_protocol_no_col: '15-111', kn.sample_id_col: 'DEV-01'},
            {kn.tm_trial_protocol_no_col: '22-222', kn.sample_id_col: 'DEV-01'},
            {kn.tm_trial_protocol_no_col: '22-222', kn.sample_id_col: 'DEV-01'}
        ]
        sort_order = sort_by_reverse_protocol_no(matches, sort_order)
        assert sort_order[('DEV-01', '11-111')] == [0, 0, 0, 0, 3]
        assert sort_order[('DEV-01', '09-999')] == [0, 0, 0, 1, 4]
        assert sort_order[('DEV-01', '15-000')] == [7, 0, 0, 0, 2]
        assert sort_order[('DEV-01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('DEV-01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('DEV-01', '15-111')] == [7, 0, 0, 0, 1]
        assert sort_order[('DEV-01', '22-222')] == [7, 1, 0, 0, 0]

    def test_final_sort(self):

        mso = {}
        sort_order = {
            ('DEV-01', '11-111'): [0, 0, 0, 0, 1],
            ('DEV-01', '09-999'): [0, 0, 0, 1, 0],
            ('DEV-01', '12-000'): [2, 3, 1, 0, 2],
            ('DEV-01', '12-222'): [2, 3, 1, 1, 3],
            ('DEV-01', '12-333'): [2, 3, 1, 0, 4],
            ('DEV-01', '13-000'): [2, 0, 0, 0, 5],
            ('DEV-01', '15-333'): [2, 3, 1, 0, 6]
        }
        mso = final_sort(sort_order, mso)
        assert mso[('DEV-01', '11-111')] == 0
        assert mso[('DEV-01', '09-999')] == 1
        assert mso[('DEV-01', '13-000')] == 2
        assert mso[('DEV-01', '12-000')] == 3
        assert mso[('DEV-01', '12-333')] == 4
        assert mso[('DEV-01', '15-333')] == 5
        assert mso[('DEV-01', '12-222')] == 6

    def test_add_sort_order(self):

        tm = get_demo_trial_matches()
        self.s = Sort(trial_matches=tm)
        tm = self.s.add_sort_order()

        # print tm[[kn.tm_trial_protocol_no_col, 'sort_order']].sort_values(by='sort_order', ascending=True)
        assert tm[[kn.tm_trial_protocol_no_col, 'sort_order']].\
            sort_values(by='sort_order', ascending=True)[kn.tm_trial_protocol_no_col].tolist() == \
            [
                '0003-000',  # tm13 (SV match (gets a sort order of -1))
                '0001-000',  # tm11 (mmr status deficient)
                '111-000',   # tm1  (tier 1, variant match, specific cancer type, DFCI, higher protocol #)
                '000-111',   # tm10 (tier 1, variant match, specific cancer type, DFCI, lower protocol #)
                '999-000',   # tm9  (tier 1, variant match, specific cancer type, MGH)
                '888-000',   # tm8  (tier 1, variant match, solid cancer type)
                '777-000',   # tm7  (tier 1, gene match)
                '444-000',   # tm4  (tier 2)
                '333-000',   # tm3  (CNV)
                '555-000',   # tm5  (tier 3)
                '666-000',   # tm6  (tier 4, higher protocol #)
                '222-000',   # tm2  (tier 4)
                '0002-000',  # tm12 (wildtype)
                '0004-000',  # tm14 (clinical only)
            ]
