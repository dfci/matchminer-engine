import unittest

from tests.test_load_service import *
from src.data_store import samples_data_model as dm
from src.services.load_service.variants_utils import VariantsUtils


class TestVariantsUtils(unittest.TestCase):

    def setUp(self):
        super(TestVariantsUtils, self).setUp()

        self.v = VariantsUtils(sample_obj={})

    def tearDown(self):
        pass

    def test_init_sample_obj(self):

        assert sorted(self.v.sample_obj.keys()) == sorted(s.gcol_list + s.signature_cols), \
            sorted(self.v.sample_obj.keys())

        for col in s.gcol_list:
            assert self.v.sample_obj[col] == [], self.v.sample_obj[col]

        for col in s.signature_cols:
            assert self.v.sample_obj[col] is None, self.v.sample_obj[col]

    def test_determine_reference_residue(self):

        ref_residue_missense = self.v._determine_reference_residue(data=mutation_missense_data)
        ref_residue_nonsense = self.v._determine_reference_residue(data=mutation_nonsense_data)

        assert ref_residue_missense == 'p.V600', ref_residue_missense
        assert ref_residue_nonsense is None, ref_residue_nonsense

        # catch missing data
        return_value_error = False
        missing_protein_change_data = mutation_missense_data.copy()
        del missing_protein_change_data[kn.protein_change_col]
        try:
            self.v._determine_reference_residue(data=missing_protein_change_data)
        except ValueError:
            return_value_error = True

        assert return_value_error is True

    def test_create_variant_object(self):

        self.v.create_variant_object(data=mutation_missense_data)
        self.v.create_variant_object(data=cnv_heterozygous_del_data)
        self.v.create_variant_object(data=sv_data)

        mutation_list = sorted(self.v.sample_obj[kn.mutation_list_col][0].keys())
        cnv_list = sorted(self.v.sample_obj[kn.cnv_list_col][0].keys())
        sv_list = sorted(self.v.sample_obj[kn.sv_list_col][0].keys())
        assert mutation_list == sorted(dm.mutations_schema.keys()), mutation_list
        assert cnv_list == sorted(dm.cnvs_schema.keys()), cnv_list
        assert sv_list == sorted(['structuralVariantComment', 'variantCategory']), sv_list

    def test_split_mmr_status(self):

        unsupported_text = 'TEST'
        res1 = self.v._split_mmr_status(None)
        res2 = self.v._split_mmr_status(s.mmr_status_cannot_assess_val)
        res3 = self.v._split_mmr_status(s.mmr_status_indeterminate_val)
        res4 = self.v._split_mmr_status('Proficient (MMR-P / MSS)')
        res5 = self.v._split_mmr_status('Deficient (MMR-D / MSI-H)')
        res6 = self.v._split_mmr_status(unsupported_text)

        assert res1 == (None, None), res1
        assert res2 == (s.mmr_status_cannot_assess_val, s.mmr_status_cannot_assess_val), res2
        assert res3 == (s.mmr_status_indeterminate_val, s.mmr_status_indeterminate_val), res3
        assert res4 == (s.mmr_status_proficient_val, s.ms_status_mss_val), res4
        assert res5 == (s.mmr_status_deficient_val, s.ms_status_msih_val), res5
        assert res6 == (unsupported_text, unsupported_text), res6

    def test_determine_signature_type(self):

        # MMR deficient
        self.v.determine_signature_type(data=signature_mmr_deficient_data)
        assert self.v.sample_obj[kn.mmr_status_col] == s.mmr_status_deficient_val, self.v.sample_obj[kn.mmr_status_col]
        assert self.v.sample_obj[kn.ms_status_col] == s.ms_status_msih_val, self.v.sample_obj[kn.ms_status_col]
        assert self.v.sample_obj[kn.apobec_status_col] is None, self.v.sample_obj[kn.apobec_status_col]

        # MMR not specified
        self.v.determine_signature_type(data=signature_mmr_none_data)
        assert self.v.sample_obj[kn.mmr_status_col] is None, self.v.sample_obj[kn.mmr_status_col]
        assert self.v.sample_obj[kn.ms_status_col] is None, self.v.sample_obj[kn.ms_status_col]

    def test_determine_wildtype(self):

        wt_gene_list = self.v.sample_obj[kn.wt_genes_col]
        self.v.determine_wildtype(data=wt1_data)
        assert sorted([i[kn.hugo_symbol_col] for i in wt_gene_list]) == sorted(['BRAF']), sorted(wt_gene_list)

        self.v.determine_wildtype(data=wt2_data)
        assert sorted([i[kn.hugo_symbol_col] for i in wt_gene_list]) == sorted(['BRAF', 'EGFR']), sorted(wt_gene_list)

        # catch missing data
        return_value_error = False
        missing_protein_change_data = {}
        try:
            self.v.determine_wildtype(data=missing_protein_change_data)
        except ValueError:
            return_value_error = True

        assert return_value_error is True

    def test_determine_low_coverage_type(self):

        pnlist = self.v.sample_obj[kn.pertinent_negatives_list_col]
        plclist = self.v.sample_obj[kn.pertinent_undercovered_list_col]
        nplclist = self.v.sample_obj[kn.additional_undercovered_list_col]

        # PN
        self.v.determine_low_coverage_type(data=pertinent_negative_data)
        assert pnlist == [pertinent_negative_data], pnlist
        assert plclist == [], plclist
        assert nplclist == [], nplclist

        # PLC
        self.v.determine_low_coverage_type(data=pertinent_undercovered_data)
        assert pnlist == [pertinent_negative_data], pnlist
        assert plclist == [pertinent_undercovered_data], plclist
        assert nplclist == [], nplclist

        # NPLC
        self.v.determine_low_coverage_type(data=additional_undercovered_data)
        assert pnlist == [pertinent_negative_data], pnlist
        assert plclist == [pertinent_undercovered_data], plclist
        assert nplclist == [additional_undercovered_data], nplclist

        # PN again
        self.v.determine_low_coverage_type(data=pertinent_negative_v2_data)
        assert sorted(pnlist) == sorted([pertinent_negative_data, pertinent_negative_v2_data]), sorted(pnlist)
        assert plclist == [pertinent_undercovered_data], plclist
        assert nplclist == [additional_undercovered_data], nplclist
