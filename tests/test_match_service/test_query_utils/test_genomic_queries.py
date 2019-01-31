from src.utilities import settings as s
from src.data_store import key_names as kn
from tests.test_match_service import TestQueryUtilitiesShared
from src.services.match_service.query_utils.genomic_queries import GenomicQueries
from src.services.match_service.query_utils.proj_utils import ProjUtils


class TestGenomicQueries(TestQueryUtilitiesShared):
    def setUp(self):
        super(TestGenomicQueries, self).setUp()

        self.gq = GenomicQueries()
        self.p = ProjUtils()
        self.standard_genomic_mut_proj = {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.cnv_call_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.sv_comment_col): 1,
            '%s.%s' % (kn.mutation_list_col, kn.tier_col): 1,
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
            '%s.%s' % (kn.cnv_list_col, kn.tier_col): 1,
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
            '%s.%s' % (kn.sv_list_col, kn.tier_col): 1,
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
            '%s.%s' % (kn.wt_genes_col, kn.tier_col): 1,
        }

    def tearDown(self):
        self.db.testSamples.drop()

    def test_create_gene_level_mutation_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_v600e,
            self.test_case_braf_non_v600e,
            self.test_case_tp53_r278w
        ])

        # create query
        q1 = self.gq.create_gene_level_query(gene_name='BRAF',
                                             variant_category=s.variant_category_mutation_val,
                                             include=True)
        eq1 = self.gq.create_gene_level_query(gene_name='BRAF',
                                              variant_category=s.variant_category_mutation_val,
                                              include=False)
        res1 = self._findalls(q1)
        eres1 = self._findalls(eq1)
        self._print(q1)
        self._print(eq1)

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_mutation_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF'])

        # assertions
        assert res1 == ['TEST-SAMPLE-BRAF-NON-V600E', 'TEST-SAMPLE-BRAF-V600E'], res1
        assert eres1 == ['TEST-SAMPLE-TP53-R278W'], eres1
        assert mr1 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF'
        }, mr1
        assert p1 == self.standard_genomic_mut_proj, p1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            kn.variant_category_col: s.variant_category_mutation_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_variant_level_mutation_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_v600e,
            self.test_case_braf_non_v600e,
            self.test_case_tp53_r278w,
            self.test_case_erbb2_v600e,
            self.test_case_no_mutation
        ])

        # BRAF V600E (inclusion)
        q1 = self.gq.create_mutation_query(gene_name='BRAF', protein_change='p.V600E', include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-V600E'], res1

        # BRAF V600E (exclusion)
        q2 = self.gq.create_mutation_query(gene_name='BRAF', protein_change='p.V600E', include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-BRAF-NON-V600E', 'TEST-SAMPLE-ERBB2-V600E',
                        'TEST-SAMPLE-NO-MUTATION', 'TEST-SAMPLE-TP53-R278W'], res2

        # BRAF V600D (inclusion)
        q3 = self.gq.create_mutation_query(gene_name='BRAF', protein_change='p.V600D', include=True)
        res3 = self._findall(q3)
        assert len(res3) == 0, res3

        # BRAF V600D (exclusion)
        q4 = self.gq.create_mutation_query(gene_name='BRAF', protein_change='p.V600D', include=False)
        res4 = self._findalls(q4)
        assert res4 == ['TEST-SAMPLE-BRAF-NON-V600E', 'TEST-SAMPLE-BRAF-V600E', 'TEST-SAMPLE-ERBB2-V600E',
                        'TEST-SAMPLE-NO-MUTATION', 'TEST-SAMPLE-TP53-R278W'], res4

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_mutation_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.protein_change_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF', 'p.V600E'])
        p3, mr3 = self.p.create_genomic_proj(include=True, query=q3)
        p4, mr4 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.protein_change_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF', 'p.V600D'])

        assert p1 == self.standard_genomic_mut_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 'p.V600E'
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.protein_change_key: 'p.V600E',
            kn.variant_category_col: s.variant_category_mutation_val
        }, p2
        assert mr2 is None, mr2
        assert p3 == self.standard_genomic_mut_proj, p3
        assert mr3 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.protein_change_col): 'p.V600D'
        }, mr3
        assert p4 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.protein_change_key: 'p.V600D',
            kn.variant_category_col: s.variant_category_mutation_val
        }, p4
        assert mr4 is None, mr4

        # clean up
        self.db.testSamples.drop()

    def test_create_wildcard_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_v600e,
            self.test_case_braf_non_v600e,
            self.test_case_tp53_r278w,
            self.test_case_erbb2_v600e,
            self.test_case_no_mutation
        ])

        # BRAF V600 wildcard (inclusion)
        q1 = self.gq.create_wildcard_query(gene_name='BRAF', protein_change='p.V600', include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-NON-V600E', 'TEST-SAMPLE-BRAF-V600E'], res1

        # BRAF V600 wildcard (exclusion)
        q2 = self.gq.create_wildcard_query(gene_name='BRAF', protein_change='p.V600', include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-ERBB2-V600E', 'TEST-SAMPLE-NO-MUTATION', 'TEST-SAMPLE-TP53-R278W'], res2

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_mutation_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.ref_residue_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF', 'p.V600'])

        assert p1 == self.standard_genomic_mut_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.ref_residue_col): 'p.V600'
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.ref_residue_key: 'p.V600',
            kn.variant_category_col: s.variant_category_mutation_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

        # clean up
        self.db.testSamples.drop()

    def test_create_variant_class_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_nonsense_mutation,
            self.test_case_braf_missense_mutation,
            self.test_case_egfr_nonsense_mutation,
            self.test_case_no_mutation
        ])

        # BRAF Nonsense (inclusion)
        q1 = self.gq.create_variant_class_query(gene_name='BRAF', variant_class='Nonsense_Mutation', include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-NONSENSE'], res1

        # BRAF Nonsense (exclusion)
        q2 = self.gq.create_variant_class_query(gene_name='BRAF', variant_class='Nonsense_Mutation', include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-BRAF-MISSENSE', 'TEST-SAMPLE-EGFR-NONSENSE', 'TEST-SAMPLE-NO-MUTATION'], res2

        # projections
        # projections
        vc = self.gq.variant_category_dict[s.variant_category_variant_class_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.variant_class_key],
                                             vals=[s.variant_category_variant_class_val, 'BRAF', 'Nonsense_Mutation'])
        assert p1 == self.standard_genomic_mut_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 'Nonsense_Mutation'
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.variant_class_key: 'Nonsense_Mutation',
            kn.variant_category_col: s.variant_category_variant_class_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_exon_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_exon_wrong_exon,
            self.test_case_exon_wrong_gene,
            self.test_case_exon_wrong_variant_class,
            self.test_case_braf_exon_20,
            self.test_case_no_mutation,
            self.test_case_erbb2_v600e
        ])

        # BRAF exon 20 (inclusion)
        q1 = self.gq.create_exon_query(gene_name='BRAF', exon=20, include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-EXON-20', 'TEST-SAMPLE-ERBB2-V600E',
                        'TEST-SAMPLE-EXON-WRONG-VARIANT-CLASS'], res1

        # BRAF exon 20 (exclusion)
        q2 = self.gq.create_exon_query(gene_name='BRAF', exon=20, include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-EXON-WRONG-EXON', 'TEST-SAMPLE-EXON-WRONG-GENE', 'TEST-SAMPLE-NO-MUTATION'], res2

        # BRAF exon 20 In Frame Insertion (inclusion)
        q3 = self.gq.create_exon_query(gene_name='BRAF', exon=20, variant_class='In_Frame_Ins', include=True)
        res3 = self._findalls(q3)
        self._print(q3)
        assert res3 == ['TEST-SAMPLE-BRAF-EXON-20'], res3

        # BRAF exon 20 In Frame Insertion (exclusion)
        q4 = self.gq.create_exon_query(gene_name='BRAF', exon=20, variant_class='In_Frame_Ins', include=False)
        res4 = self._findalls(q4)
        self._print(q4)
        assert res4 == ['TEST-SAMPLE-ERBB2-V600E', 'TEST-SAMPLE-EXON-WRONG-EXON', 'TEST-SAMPLE-EXON-WRONG-GENE',
                        'TEST-SAMPLE-EXON-WRONG-VARIANT-CLASS', 'TEST-SAMPLE-NO-MUTATION'], res4

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_mutation_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.transcript_exon_key, self.p.variant_class_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF', 20, None])
        p3, mr3 = self.p.create_genomic_proj(include=True, query=q3)
        p4, mr4 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.transcript_exon_key, self.p.variant_class_key],
                                             vals=[s.variant_category_mutation_val, 'BRAF', 20, 'In_Frame_Ins'])

        assert p1 == self.standard_genomic_mut_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 20
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.transcript_exon_key: 20,
            kn.variant_category_col: s.variant_category_mutation_val
        }, p2
        assert mr2 is None, mr2
        assert p3 == self.standard_genomic_mut_proj, p3
        assert mr3 == {
            '%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.mutation_list_col, kn.transcript_exon_col): 20,
            '%s.%s' % (kn.mutation_list_col, kn.variant_class_col): 'In_Frame_Ins'
        }, mr3
        assert p4 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.transcript_exon_key: 20,
            self.p.variant_class_key: 'In_Frame_Ins',
            kn.variant_category_col: s.variant_category_mutation_val
        }, p4
        assert mr4 is None, mr4

        # clean up
        self.db.testSamples.drop()

    def test_create_gene_level_cnv_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_generic_cnv,
            self.test_case_braf_cnv_hetero_del,
            self.test_case_braf_cnv_gain,
            self.test_case_no_cnv,
            self.test_case_braf_v600e
        ])

        # BRAF any CNV (inclusion)
        q1 = self.gq.create_gene_level_query(gene_name='BRAF',
                                             variant_category=s.variant_category_cnv_val,
                                             include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-CNV-GAIN', 'TEST-SAMPLE-BRAF-CNV-HETERO-DEL',
                        'TEST-SAMPLE-BRAF-GENERIC-CNV'], res1

        # BRAF any CNV (exclusion)
        q2 = self.gq.create_gene_level_query(gene_name='BRAF',
                                             variant_category=s.variant_category_cnv_val,
                                             include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-BRAF-V600E', 'TEST-SAMPLE-NO-CNV'], res2

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_cnv_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key],
                                             vals=[s.variant_category_cnv_val, 'BRAF'])
        assert p1 == self.standard_genomic_cnv_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF',
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            kn.variant_category_col: s.variant_category_cnv_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_cnv_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_cnv_hetero_del,
            self.test_case_braf_cnv_gain,
            self.test_case_no_cnv,
            self.test_case_braf_v600e
        ])

        # BRAF CNV Heterozygous deletion (inclusion)
        q1 = self.gq.create_cnv_query(gene_name='BRAF', cnv_call=s.cnv_call_hetero_del, include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-CNV-HETERO-DEL'], res1

        # BRAF CNV Heterozygous deletion (exclusion)
        q2 = self.gq.create_cnv_query(gene_name='BRAF', cnv_call=s.cnv_call_hetero_del, include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-BRAF-CNV-GAIN', 'TEST-SAMPLE-BRAF-V600E', 'TEST-SAMPLE-NO-CNV'], res2

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_cnv_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key, self.p.cnv_call_key],
                                             vals=[s.variant_category_cnv_val, 'BRAF', s.cnv_call_hetero_del])
        assert p1 == self.standard_genomic_cnv_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF',
            '%s.%s' % (kn.cnv_list_col, kn.cnv_call_col): s.cnv_call_hetero_del,
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            self.p.cnv_call_key: s.cnv_call_hetero_del,
            kn.variant_category_col: s.variant_category_cnv_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_sv_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_sv,
            self.test_case_sv_2,
            self.test_case_no_sv
        ])

        # BRAF SV (inclusion)
        q1 = self.gq.create_sv_query(gene_name='NTRK1', include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-NTRK1-SV'], res1

        # BRAF SV (exclusion)
        q2 = self.gq.create_sv_query(gene_name='NTRK1', include=False)
        res2 = self._findalls(q2)
        self._print(q2)
        assert res2 == ['TEST-SAMPLE-NO-SV', 'TEST-SAMPLE-NTRK2-SV'], res2

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_sv_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key],
                                             vals=[s.variant_category_sv_val, 'NTRK1'])
        assert p1 == self.standard_genomic_sv_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.sv_list_col, kn.sv_comment_col): q1[kn.sv_list_col]['$elemMatch'][kn.sv_comment_col]
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'NTRK1',
            kn.variant_category_col: s.variant_category_sv_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_gene_level_wt_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_wt
        ])

        # BRAF WT (inclusion)
        q1 = self.gq.create_gene_level_query(gene_name='BRAF',
                                             variant_category=s.variant_category_wt_val,
                                             include=True)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-BRAF-WT']

        # BRAF ST (exclusion)
        q2 = self.gq.create_gene_level_query(gene_name='BRAF',
                                             variant_category=s.variant_category_wt_val,
                                             include=False)
        res2 = self._findall(q2)
        self._print(q2)
        assert res2 == []

        # projections
        vc = self.gq.variant_category_dict[s.variant_category_wt_val]
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        p2, mr2 = self.p.create_genomic_proj(include=False,
                                             keys=[vc, self.p.hugo_symbol_key],
                                             vals=[s.variant_category_wt_val, 'BRAF'])
        assert p1 == self.standard_genomic_wt_proj, p1
        assert mr1 == {
            '%s.%s' % (kn.wt_genes_col, kn.hugo_symbol_col): 'BRAF'
        }, mr1
        assert p2 == {
            self.p.hugo_symbol_key: 'BRAF',
            kn.variant_category_col: s.variant_category_wt_val
        }, p2
        assert mr2 is None, mr2

        # clean up
        self.db.testSamples.drop()

    def test_create_mutational_signature_query(self):

        # set up
        self.db.testSamples.insert(self.test_case_mmr_deficient)

        # MMR Deficient (inclusion)
        q1 = self.gq.create_mutational_signature_query(signature_type=kn.mmr_status_col,
                                                       signature_val=s.mmr_status_deficient_val)
        res1 = self._findalls(q1)
        self._print(q1)
        assert res1 == ['TEST-SAMPLE-MMR-DEFICIENT'], res1

        # MMR Deficient (exclusion)
        q2 = self.gq.create_mutational_signature_query(signature_type=kn.mmr_status_col,
                                                       signature_val=s.mmr_status_proficient_val)
        res2 = self._findall(q2)
        self._print(q2)
        assert len(res2) == 0, res2

        # projections
        p1, mr1 = self.p.create_genomic_proj(include=True, query=q1)
        assert p1 == {
            '_id': 0, kn.sample_id_col: 1, kn.mrn_col: 1, kn.vital_status_col: 1,
            kn.mmr_status_col: 1
        }, p1
        assert mr1 == {
            kn.mmr_status_col: s.mmr_status_deficient_val
        }

        # clean up
        self.db.testSamples.drop()

    def test_create_any_variant_query(self):

        # set up
        self.db.testSamples.insert_many([
            self.test_case_braf_v600e,
            self.test_case_erbb2_v600e,
            self.test_case_tp53_r278w,
            self.test_case_no_mutation,
            self.test_case_braf_generic_cnv
        ])

        # todo deprecated
        # # BRAF Any Variation (inclusion)
        # q1 = self.gq.create_any_variant_query(gene_name='BRAF', include=True)
        # res1 = self._findalls(q1)
        # self._print(q1)
        # assert res1 == ['TEST-SAMPLE-BRAF-GENERIC-CNV', 'TEST-SAMPLE-BRAF-V600E', 'TEST-SAMPLE-ERBB2-V600E']
        #
        # # BRAF Any Variation (exclusion)
        # q2 = self.gq.create_any_variant_query(gene_name='BRAF', include=False)
        # res2 = self._findalls(q2)
        # self._print(q2)
        # assert res2 == ['TEST-SAMPLE-NO-MUTATION', 'TEST-SAMPLE-TP53-R278W']
        #
        # # projections
        # p1, mr1 = self.p.create_any_variant_proj(gene_name='BRAF', include=True, query=q1)
        # p2, mr2 = self.p.create_any_variant_proj(gene_name='BRAF', include=False, query=q1)
        # assert p1 == {
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
        # }, p1
        # assert mr1 == {
        #     '$or': [
        #         {'%s.%s' % (kn.mutation_list_col, kn.hugo_symbol_col): 'BRAF'},
        #         {'%s.%s' % (kn.cnv_list_col, kn.hugo_symbol_col): 'BRAF'}
        #     ]
        # }, mr1
        # self._print(p2)
        # assert p2 == {
        #     '$and': [
        #         {self.p.hugo_symbol_key: 'BRAF', kn.variant_category_col: s.variant_category_mutation_val},
        #         {self.p.hugo_symbol_key: 'BRAF', kn.variant_category_col: s.variant_category_cnv_val}
        #     ]
        # }, p2
        # assert mr2 is None, mr2
        #
        # # clean up
        # self.db.testSamples.drop()

    def test_create_low_coverage_query(self):
        pass
