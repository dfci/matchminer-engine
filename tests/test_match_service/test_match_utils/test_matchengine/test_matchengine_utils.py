from tests.test_match_service import TestQueryUtilitiesShared
from src.services.match_service.match_utils.matchengine.matchengine_utils import *


class TestMatchEngineUtils(TestQueryUtilitiesShared):

    def setUp(self):
        super(TestMatchEngineUtils, self).setUp()

    def tearDown(self):
        pass

    def test_assess_inclusion(self):
        assert assess_inclusion('Mutation') is True
        assert assess_inclusion('!Mutation') is False

    def test_sanitize_exclusion_vals(self):
        assert sanitize_exclusion_vals('!Mutation') == 'Mutation'
        assert sanitize_exclusion_vals('Mutation') == 'Mutation'

    def test_normalize_variant_category_val(self):
        assert normalize_variant_category_val(s.mt_cnv_val) == s.variant_category_cnv_val
        assert normalize_variant_category_val(s.mt_mut_val) == s.variant_category_mutation_val
        assert normalize_variant_category_val(s.mt_sv_val) == s.variant_category_sv_val
        assert normalize_variant_category_val(s.mt_any_vc_val) == s.variant_category_any_val

    def test_normalize_cnv_call_val(self):
        assert normalize_cnv_call_val(s.mt_high_amp_val) == s.cnv_call_high_amp
        assert normalize_cnv_call_val(s.mt_homo_del_val) == s.cnv_call_homo_del
        assert normalize_cnv_call_val(s.mt_hetero_del_val) == s.cnv_call_hetero_del
        assert normalize_cnv_call_val(s.mt_low_amp_val) == s.cnv_call_gain

    def test_normalize_signature_vals(self):
        assert normalize_signature_vals(s.mt_mmr_status, s.mt_mmr_deficient_val) == (kn.mmr_status_col, s.mmr_status_deficient_val)
        assert normalize_signature_vals(s.mt_mmr_status, s.mt_mmr_proficient_val) == (kn.mmr_status_col, s.mmr_status_proficient_val)
        assert normalize_signature_vals(s.mt_ms_status, s.mt_msi_high_val) == (kn.ms_status_col, s.ms_status_msih_val)
        assert normalize_signature_vals(s.mt_ms_status, s.mt_mss_val) == (kn.ms_status_col, s.ms_status_mss_val)
        assert normalize_signature_vals(s.mt_tobacco_status, 'Yes') == (kn.tobacco_status_col, 'Yes')
        assert normalize_signature_vals(s.mt_tmz_status, 'Yes') == (kn.tmz_status_col, 'Yes')
        assert normalize_signature_vals(s.mt_pole_status, 'Yes') == (kn.pole_status_col, 'Yes')
        assert normalize_signature_vals(s.mt_apobec_status, 'Yes') == (kn.apobec_status_col, 'Yes')
        assert normalize_signature_vals(s.mt_uva_status, 'Yes') == (kn.uva_status_col, 'Yes')

    def test_reformat_matches(self):
        m1 = [
            {kn.mutation_list_col: {kn.hugo_symbol_col: 'BRAF'}},
            {kn.mutation_list_col: {kn.hugo_symbol_col: 'EGFR'}},
            {kn.cnv_list_col: {kn.hugo_symbol_col: 'AKT1'}},
            {kn.wt_genes_col: {kn.hugo_symbol_col: 'TP53'}},
            {kn.sv_list_col: {kn.sv_comment_col: 'NTRK1'}}
        ]
        assert sorted(reformat_matches(m1)) == sorted([
            {kn.mutation_list_col: [{kn.hugo_symbol_col: 'BRAF'}]},
            {kn.mutation_list_col: [{kn.hugo_symbol_col: 'EGFR'}]},
            {kn.cnv_list_col: [{kn.hugo_symbol_col: 'AKT1'}]},
            {kn.wt_genes_col: [{kn.hugo_symbol_col: 'TP53'}]},
            {kn.sv_list_col: [{kn.sv_comment_col: 'NTRK1'}]}
        ])

