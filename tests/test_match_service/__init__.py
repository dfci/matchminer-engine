import os
import json
import yaml
import unittest
import pandas as pd
import datetime as dt

from src.utilities import settings as s
s.MONGO_URI = 'mongodb://localhost:27017'
s.MONGO_DBNAME = 'matchminer'

from src.utilities.utilities import get_db
from src.data_store import key_names as kn
from src.data_store.samples_data_model import vital_status_allowed_vals
from src.data_store.trial_matches_data_model import trial_matches_schema


class TestQueryUtilitiesShared(unittest.TestCase):

    def setUp(self):

        self.db = get_db(mongo_uri=s.MONGO_URI, mongo_dbname=s.MONGO_DBNAME)
        self.data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../', 'data'))
        self.proj = {
            '_id': 0,
            kn.sample_id_col: 1,
            kn.mrn_col: 1,
            kn.vital_status_col: 1
        }

        # test data
        self.test_case_lung = {
            kn.sample_id_col: 'TEST-SAMPLE-LUNG',
            kn.oncotree_primary_diagnosis_name_col: 'Lung'
        }
        self.test_case_colon = {
            kn.sample_id_col: 'TEST-SAMPLE-COLON',
            kn.oncotree_primary_diagnosis_name_col: 'Colon'
        }
        self.test_case_child = {
            kn.sample_id_col: 'TEST-SAMPLE-CHILD',
            kn.birth_date_col: dt.datetime.today() - dt.timedelta(days=360)
        }
        self.test_case_adult = {
            kn.sample_id_col: 'TEST-SAMPLE-ADULT',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1)
        }
        self.test_case_male = {
            kn.sample_id_col: 'TEST-SAMPLE-MALE',
            kn.gender_col: 'Male'
        }
        self.test_case_female = {
            kn.sample_id_col: 'TEST-SAMPLE-FEMALE',
            kn.gender_col: 'Female'
        }
        self.test_case_braf_v600e = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-V600E',
            kn.oncotree_primary_diagnosis_name_col: 'Breast',
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.protein_change_col: 'p.V600E',
                    kn.ref_residue_col: 'p.V600'
                },
                {
                    kn.hugo_symbol_col: 'ERBB2',
                    kn.protein_change_col: 'p.D769H',
                    kn.ref_residue_col: 'p.D769'
                }
            ]
        }
        self.test_case_erbb2_v600e = {
            kn.sample_id_col: 'TEST-SAMPLE-ERBB2-V600E',
            kn.mutation_list_col: [
                {
                    kn.hugo_symbol_col: 'ERBB2',
                    kn.protein_change_col: 'p.V600E',
                    kn.ref_residue_col: 'p.V600',
                    kn.transcript_exon_col: 19,
                },
                {
                    kn.hugo_symbol_col: 'BRAF',
                    kn.protein_change_col: 'p.D769H',
                    kn.ref_residue_col: 'p.D769',
                    kn.transcript_exon_col: 20,
                }
            ]
        }
        self.test_case_tp53_r278w = {
            kn.sample_id_col: 'TEST-SAMPLE-TP53-R278W',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'TP53',
                kn.protein_change_col: 'p.R278W',
                kn.ref_residue_col: 'pp.R278'
            }]
        }
        self.test_case_braf_non_v600e = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-NON-V600E',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600X',
                kn.ref_residue_col: 'p.V600'
            }]
        }
        self.test_case_no_mutation = {
            kn.sample_id_col: 'TEST-SAMPLE-NO-MUTATION',
            kn.mutation_list_col: []
        }
        self.test_case_braf_exon_20 = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-EXON-20',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.transcript_exon_col: 20,
                kn.variant_class_col: 'In_Frame_Ins'
            }]
        }
        self.test_case_exon_wrong_variant_class = {
            kn.sample_id_col: 'TEST-SAMPLE-EXON-WRONG-VARIANT-CLASS',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.transcript_exon_col: 20,
                kn.variant_class_col: 'In_Frame_Del'
            }]
        }
        self.test_case_exon_wrong_exon = {
            kn.sample_id_col: 'TEST-SAMPLE-EXON-WRONG-EXON',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.transcript_exon_col: 19,
                kn.variant_class_col: 'In_Frame_Ins'
            }]
        }
        self.test_case_exon_wrong_gene = {
            kn.sample_id_col: 'TEST-SAMPLE-EXON-WRONG-GENE',
            kn.mutation_list_col: [{
                kn.hugo_symbol_col: 'EGFR',
                kn.transcript_exon_col: 20,
                kn.variant_class_col: 'In_Frame_Ins'
            }]
        }
        self.test_case_braf_generic_cnv = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-GENERIC-CNV',
            kn.oncotree_primary_diagnosis_name_col: 'Hodgkin Lymphoma',
            kn.birth_date_col: dt.datetime(year=1900, day=1, month=1),
            kn.cnv_list_col: [{
                kn.hugo_symbol_col: 'BRAF'
            }]
        }
        self.test_case_braf_cnv_hetero_del = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-CNV-HETERO-DEL',
            kn.oncotree_primary_diagnosis_name_col: 'Hodgkin Lymphoma',
            kn.birth_date_col: dt.datetime.today() - dt.timedelta(days=360),
            kn.cnv_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.cnv_call_col: s.cnv_call_hetero_del
            }]
        }
        self.test_case_braf_cnv_gain = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-CNV-GAIN',
            kn.cnv_list_col: [{
                kn.hugo_symbol_col: 'BRAF',
                kn.cnv_call_col: s.cnv_call_gain
            }]
        }
        self.test_case_no_cnv = {
            kn.sample_id_col: 'TEST-SAMPLE-NO-CNV',
            kn.cnv_list_col: []
        }
        self.test_case_sv = {
            kn.sample_id_col: 'TEST-SAMPLE-NTRK1-SV',
            kn.sv_list_col: [{
                kn.sv_comment_col: 'This sample has a NTRK1 fusion'
            }]
        }
        self.test_case_sv_2 = {
            kn.sample_id_col: 'TEST-SAMPLE-NTRK2-SV',
            kn.sv_list_col: [{
                kn.sv_comment_col: 'This sample has a NTRK2 fusion'
            }]
        }
        self.test_case_no_sv = {
            kn.sample_id_col: 'TEST-SAMPLE-NO-SV',
            kn.sv_list_col: []
        }
        self.test_case_braf_wt = {
            kn.sample_id_col: 'TEST-SAMPLE-BRAF-WT',
            kn.wt_genes_col: [{
                kn.hugo_symbol_col: 'BRAF'
            }]
        }

        self.test_case_mmr_deficient = {
            kn.sample_id_col: 'TEST-SAMPLE-MMR-DEFICIENT',
            kn.mmr_status_col: s.mmr_status_deficient_val
        }

        self.test_cases = [
            self.test_case_lung,
            self.test_case_colon,
            self.test_case_braf_v600e,
            self.test_case_braf_non_v600e,
            self.test_case_no_mutation,
            self.test_case_braf_generic_cnv,
            self.test_case_braf_cnv_hetero_del,
            self.test_case_braf_cnv_gain,
            self.test_case_sv,
            self.test_case_sv_2,
            self.test_case_mmr_deficient,
            self.test_case_braf_wt,
            self.test_case_braf_exon_20
        ]

        # test values
        self.all_nsclc_cancer_types = sorted([
            'Adenoid Cystic Carcinoma of the Lung',
            'Giant Cell Carcinoma of the Lung',
            'Large Cell Lung Carcinoma With Rhabdoid Phenotype',
            'Poorly Differentiated Non-Small Cell Lung Cancer',
            'Lung Adenosquamous Carcinoma',
            'Basaloid Large Cell Carcinoma of the Lung',
            'Large Cell Lung Carcinoma',
            'Spindle Cell Carcinoma of the Lung',
            'Non-Small Cell Lung Cancer',
            'Lymphoepithelioma-like Carcinoma of the Lung',
            'Mucoepidermoid Carcinoma of the Lung',
            'Lung Adenocarcinoma',
            'Lung Squamous Cell Carcinoma',
            'Salivary Gland-Type Tumor of the Lung',
            'Clear Cell Carcinoma of the Lung'
        ])
        self.all_liquid_cancer_types = sorted([
            'Angioimmunoblastic T-Cell Lymphoma',
            'Acute Lymphoid Leukemia',
            'Mycosis Fungoides',
            'Diffuse Large B-Cell Lymphoma',
            'Splenic Marginal Zone Lymphoma',
            'Acute Myeloid Leukemia',
            'Chronic Myelogenous Leukemia',
            'Burkitt Lymphoma',
            'Waldenstrom Macroglobulinemia',
            'Chronic Lymphocytic Leukemia',
            'Extranodal Marginal Zone Lymphoma',
            'Peripheral T-Cell Lymphoma',
            'Systemic Mastocytosis',
            'Mantle Cell Lymphoma',
            'Peripheral T-Cell Lymphoma, NOS',
            'Myelodysplasia',
            'Hairy Cell Leukemia',
            'Marginal Zone Lymphoma',
            'Non-Langerhans Cell Histiocytosis/Erdheim-Chester Disease',
            'Hodgkin Lymphoma',
            'Chronic Myelomonocytic Leukemia',
            'Small Lymphocytic Lymphoma',
            'Leukemia',
            'B-Cell Acute Lymphoid Leukemia',
            'Lymph',
            'T-Cell Acute Lymphoid Leukemia',
            'Nodular Lymphocyte-Predominant Hodgkin Lymphoma',
            'Mediastinal Large B-Cell Lymphoma',
            'Blood',
            'Nodal Marginal Zone Lymphoma',
            'Sezary Syndrome',
            'Non-Hodgkin Lymphoma',
            'Rosai-Dorfman Disease',
            'Multiple Myeloma',
            'Polycythemia Vera',
            'T-Cell and Natural Killer Lymphoma',
            'Langerhans Cell Histiocytosis',
            'Primary CNS Lymphoma',
            'Anaplastic Large Cell Lymphoma',
            'Large Granular Lymphocytic Leukemia',
            'Follicular Lymphoma',
            'Classical Hodgkin Lymphoma',
            'Blastic Plasmacytoid Dendritic Cell Neoplasm',
            'Primary Effusion Lymphoma',
            'Essential Thrombocythaemia',
            'Acute Monocytic Leukemia',
            'B-Cell Lymphoma',
            'Myeloproliferative Neoplasm',
            'Cutaneous T-Cell Lymphoma',
            'Myelofibrosis/Osteomyelofibrosis',
            'Histiocytosis'])

        self.all_solid_cancer_types = sorted([
            'Parosteal Osteosarcoma',
            'Esophageal Adenocarcinoma',
            'Upper Tract Urothelial Carcinoma',
            'Endometrioid Borderlin Ovarian Tumor',
            'Spindle Cell Rhabdomyosarcoma',
            'Uterine Endometrioid Carcinoma',
            'Atypical Meningioma',
            'Clear Cell Sarcoma',
            'Mucosal Melanoma of the Esophagus',
            'Phyllodes Tumor of the Breast',
            'Melanotic Medulloblastoma',
            'Breast Angiosarcoma',
            'Small Cell Carcinoma of Unknown Primary',
            'Sertoli-Leydig Cell Tumor',
            'Atypical Lung Carcinoid',
            'Urethral Adenocarcinoma',
            'Combined Small Cell Lung Carcinoma',
            'Pineocytoma',
            'Mesenchymal Chondrosarcoma of the CNS',
            'Breast Invasive Lobular Carcinoma',
            'Nasopharyngeal Carcinoma',
            'Sebaceous Carcinoma',
            'Lymphoepithelioma-like Carcinoma of the Lung',
            'Conventional Type Chordoma',
            'Urethral Urothelial Carcinoma',
            'Low-Grade Fibromyxoid Sarcoma',
            'Anaplastic Astrocytoma',
            'Odontogenic Carcinoma',
            'Solid Papillary Carcinoma of the Breast',
            'Well-Differentiated Neuroendocrine Tumor of the Appendix',
            'Verrucous Penile Squamous Cell Carcinoma',
            'Acinic Cell Carcinoma',
            'Chordoma',
            'Central Neurocytoma',
            'Sellar Tumor',
            'Malignant Tumor',
            'Sarcomatoid Carcinoma of the Lung',
            'Uterine Adenosquamous Carcinoma',
            'Intestinal Type Stomach Adenocarcinoma',
            'Ossifying Fibromyxoid Tumor',
            'Pineal Tumor',
            'Astrocytoma',
            'Low-Grade Neuroepithelial Tumor',
            'Microcystic Adnexal Carcinoma',
            'Renal Mucinous Tubular Spindle Cell Carcinoma',
            'Papillary Glioneuronal Tumor',
            'Stomach Adenocarcinoma',
            'High-Grade Serous Ovarian Cancer',
            'Poorly Differentiated Carcinoma, NOS',
            'Pleomorphic Xanthoastrocytoma',
            'Gangliocytoma',
            'Esophagus/Stomach',
            'Poorly Differentiated Carcinoma of the Stomach',
            'NUT Midline Carcinoma of the Head and Neck',
            'No OncoTree Node Found',
            'Small Cell Lung Cancer',
            'Adenocarcinoma, NOS',
            'Prostate Neuroendocrine Carcinoma',
            'Radiation-Associated Sarcoma',
            'Ovarian Germ Cell Tumor',
            'Dermatofibroma',
            'Uterine Leiomyosarcoma',
            'Renal Cell Carcinoma',
            'Lentigo Maligna Melanoma',
            'Squamous Cell Carcinoma of the Vulva/Vagina',
            'Low-Grade Glioma, NOS',
            'Endocervical Adenocarcinoma',
            'Plasmacytoid/Signet Ring Cell Bladder Carcinoma',
            'Chromophobe Renal Cell Carcinoma',
            'Thymus',
            'Yolk Sac Tumor',
            'Extramammary Paget Disease',
            'Small Cell Osteosarcoma',
            'Inflammatory Myofibroblastic Bladder Tumor',
            'Atypical Teratoid/Rhabdoid Tumor',
            'Conjunctival Melanoma',
            'Spiroma/Spiradenoma',
            'Astroblastoma',
            'Seminoma',
            'Complete Hydatidiform Mole',
            'Desmoplastic Infantile Astrocytoma',
            'Dedifferentiated Liposarcoma',
            'Head and Neck',
            'Ependymoma',
            'Metaplastic Squamous Cell Carcinoma',
            'Oropharynx Squamous Cell Carcinoma',
            'Melanoma',
            'Nerve Sheath Tumor',
            'Breast Invasive Mixed Mucinous Carcinoma',
            'Testicular Lymphoma',
            'Benign Phyllodes Tumor of the Breast',
            'Head and Neck Squamous Cell Carcinoma',
            'Melanoma of Unknown Primary',
            'Granular Cell Tumor',
            'Renal Medullary Carcinoma',
            'Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma',
            'Primitive Neuroectodermal Tumor',
            'Pheochromocytoma',
            'Thymic Epithelial Tumor',
            'Ocular Melanoma',
            'Hemangiopericytoma of the Central Nervous System',
            'Ependymomal Tumor',
            'Choriocarcinoma',
            'Liver',
            'CNS/Brain',
            'Myxoid/Round-Cell Liposarcoma',
            'Myxopapillary Ependymoma',
            'Head and Neck Squamous Cell Carcinoma of Unknown Primary',
            'Dysembryoplastic Neuroepithelial Tumor',
            'Renal Clear Cell Carcinoma',
            'Gastric Type Mucinous Carcinoma',
            'Glioblastoma Multiforme',
            'Chondroblastic Osteosarcoma',
            'Miscellaneous Brain Tumor',
            'Ovarian Seromucinous Carcinoma',
            'Urethral Squamous Cell Carcinoma',
            'Prostate Sqamous Cell Carcinoma',
            'Cholangiocarcinoma',
            'Endocrine Mucin Producing Sweat Gland Carcinoma',
            'Desmoid/Aggressive Fibromatosis',
            'Lung',
            'Atypical Fibroxanthoma',
            'Malignant Peripheral Nerve Sheath Tumor',
            'Lung Squamous Cell Carcinoma',
            'root',
            'Fibroblastic Osteosarcoma',
            'Gallbladder Cancer',
            'Porphyria Cutania Tarda',
            'Duodenal Adenocarcinoma',
            'Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor',
            'Soft Tissue Myoepithelial Carcinoma',
            'Head and Neck Carcinoma, Other',
            'Liver Angiosarcoma',
            'Malignant Teratoma',
            'Vulva/Vagina',
            'Intraductal Papillary Mucinous Neoplasm',
            'Solid Pseudopapillary Neoplasm of the Pancreas',
            'Uterine Perivascular Epithelioid Cell Tumor',
            'Oligodendroglioma',
            'Undifferentiated Carcinoma of the Pancreas',
            'Pleural Mesothelioma, Epithelioid Type',
            'Invasive Breast Carcinoma',
            'Anorectal Mucosal Melanoma',
            'Basaloid Large Cell Carcinoma of the Lung',
            'Sinonasal Squamous Cell Carcinoma',
            'Renal Angiomyolipoma',
            'Prostate Adenocarcinoma',
            'Primary Neuroepithelial Tumor',
            'Breast Lobular Carcinoma In Situ',
            'Fibrolamellar Carcinoma',
            'Papillary Renal Cell Carcinoma',
            'Choroid Plexus Carcinoma',
            'Myxofibrosarcoma',
            'Non-Small Cell Lung Cancer',
            'Ovarian Epithelial Tumor',
            'Cellular Schwannoma',
            'Clear Cell Papillary Renal Cell Carcinoma',
            'Sweat Gland Carcinoma/Apocrine Eccrine Carcinoma',
            'Oral Cavity Squamous Cell Carcinoma',
            'Eye',
            'Follicular Dendritic Cell Sarcoma',
            'Gastric Remnant Adenocarcinoma',
            'Esophageal Poorly Differentiated Carcinoma',
            'Cancer of Unknown Primary',
            'Solitary Fibrous Tumor/Hemangiopericytoma',
            'Neuroendocrine Tumor, NOS',
            'Bladder/Urinary Tract',
            'Clear Cell Borderline Ovarian Tumor',
            'Breast Sarcoma',
            'Cutaneous Squamous Cell Carcinoma',
            'Embryonal Tumor',
            'Sclerosing Epithelioid Fibrosarcoma',
            'Breast Invasive Cancer, NOS',
            'Urachal Adenocarcinoma',
            'Adenomyoepithelioma of the Breast',
            'Subependymoma',
            'Periosteal Osteosarcoma',
            'Low-Grade Central Osteosarcoma',
            'Carcinoma with Osseous Metaplasia',
            'Well-Differentiated Thyroid Cancer',
            'Endometrioid Ovarian Cancer',
            'Breast Mixed Ductal and Lobular Carcinoma',
            'Germ Cell Tumor of the Vulva',
            'Cutaneous Melanoma',
            'Histiocytic Dendritic Cell Sarcoma',
            'Spindle Cell Oncocytoma of the Adenohypophysis',
            'Liposarcoma',
            'Desmoplastic/Nodular Medulloblastoma',
            'Myofibroma',
            'Small Cell Carcinoma of the Ovary',
            'Colorectal Adenocarcinoma',
            'Cervix',
            'Testis',
            'Thymoma',
            'Invasive Hydatidiform Mole',
            'Villoglandular Carcinoma',
            'Cervical Leiomyosarcoma',
            'Pineoblastoma',
            'Uterine Undifferentiated Carcinoma',
            'Poorly Differentiated Carcinoma of the Uterus',
            'Mammary Analogue Secretory Carcinoma of Salivary Gland Origin',
            'Metaplastic Carcinosarcoma',
            'Schwannoma',
            'Undifferentiated Malignant Neoplasm',
            'Rectal Adenocarcinoma',
            'Breast',
            'Uterine Sarcoma/Mesenchymal',
            'Anaplastic Ependymoma',
            'Breast Invasive Ductal Carcinoma',
            'Paraganglioma',
            'Well-Differentiated Neuroendocrine Tumor of the Rectum',
            'Medullary Carcinoma of the Colon',
            'Solitary Fibrous Tumor of the Central Nervous System',
            'Desmoplastic Trichoepithelioma',
            'Atypical Pituitary Adenoma',
            'Translocation-Associated Renal Cell Carcinoma',
            'Vaginal Adenocarcinoma',
            'Anaplastic Oligoastrocytoma',
            'Acral Melanoma',
            'Bladder Urothelial Carcinoma',
            'Small Bowel Cancer',
            'Teratoma with Malignant Transformation',
            'Ovarian Seromucinous Borderline Tumor',
            'Spindle Cell/Sclerosing Rhabdomyosarcoma',
            'Partial Hydatidiform Mole',
            'Leiomyosarcoma',
            'Uterine Dedifferentiated Carcinoma',
            'Medullary Thyroid Cancer',
            'Collecting Duct Renal Cell Carcinoma',
            'Cerebellar Liponeurocytoma',
            'Rhabdoid Cancer',
            'Medulloepithelioma',
            'Gonadoblastoma',
            'Malignant Phyllodes Tumor of the Breast',
            'Small Cell Carcinoma of the Cervix',
            'Large Cell Lung Carcinoma With Rhabdoid Phenotype',
            'Pseudomyogenic Hemangioendothelioma',
            'Cervical Adenocarcinoma',
            'Head and Neck Neuroendocrine Carcinoma',
            'Hyalinizing Trabecular Adenoma of the Thyroid',
            'Desmoplastic Melanoma',
            'Mucinous Stomach Adenocarcinoma',
            'Pancreas',
            'Clear Cell Odontogenic Carcinoma',
            'Warty Penile Squamous Cell Carcinoma',
            'Craniopharyngioma, Adamantinomatous Type',
            'Signet Ring Cell Carcinoma of the Stomach',
            'Small Cell Glioblastoma',
            'Primary Brain Tumor',
            'Cervical Adenosquamous Carcinoma',
            'Uterus',
            'Uterine Smooth Muscle Tumor',
            'Mature Teratoma',
            'Melanotic Schwannoma',
            'Serous Borderline Ovarian Tumor, Micropapillary',
            'Uterine Adenosarcoma',
            'Embryonal Rhabdomyosarcoma',
            'Giant Cell Tumor of Bone',
            'Papillary Stomach Adenocarcinoma',
            'Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma',
            'Paget Disease of the Nipple',
            'Uterine Mixed Endometrial Carcinoma',
            'Squamous Cell Carcinoma, NOS',
            'Pituitary Carcinoma',
            'Clear Cell Carcinoma of the Lung',
            'Pituicytoma',
            'Epithelioid Trophoblastic Tumor',
            'Acinar Cell Carcinoma of the Pancreas',
            'Sex Cord Stromal Tumor',
            'Mucinous Carcinoma',
            'Other',
            'Dedifferentiated Chordoma',
            'Angiosarcoma',
            'Molar Pregnancy',
            'Sarcomatoid Carcinoma of the Urinary Bladder',
            'Cervical Squamous Cell Carcinoma',
            'Diffuse Intrinsic Pontine Glioma',
            'Fibrosarcoma',
            'Epithelioid Sarcoma',
            'Extrahepatic Cholangiocarcinoma',
            'Large Cell Lung Carcinoma',
            'Teratoma',
            'Ovary/Fallopian Tube',
            'Encapsulated Glioma',
            'Anaplastic Oligodendroglioma',
            'Carcinoma with Chondroid Metaplasia',
            'Choroid Plexus Tumor',
            'Goblet Cell Carcinoid of the Appendix',
            'Placental Site Trophoblastic Tumor',
            'Dermatofibrosarcoma Protuberans',
            'Uterine Leiomyoma',
            'Inflammatory Breast Cancer',
            'Alveolar Soft Part Sarcoma',
            'Meningothelial Tumor',
            'Embryonal Carcinoma',
            'Bone',
            'Chordoid Glioma of the Third Ventricle',
            'Desmoplastic Small-Round-Cell Tumor',
            'Salivary Carcinoma',
            'Adrenal Gland',
            'Signet Ring Cell Adenocarcinoma of the Colon and Rectum',
            'Extraventricular Neurocytoma',
            'Small Bowel Well-Differentiated Neuroendocrine Tumor',
            'High-Grade Endometrial Stromal Sarcoma',
            'Hemangioblastoma',
            'Acinar Cell Carcinoma, NOS',
            'Mesonephric Carcinoma',
            'Mucinous Cystic Neoplasm',
            'Mucinous Adenocarcinoma of the Colon and Rectum',
            'Cystic Tumor of the Pancreas',
            'Bladder Adenocarcinoma',
            'Synovial Sarcoma',
            'Endometrial Carcinoma',
            'Osteoblastic Osteosarcoma',
            'Brenner Tumor, Benign',
            'Malignant Nonepithelial Tumor of the Liver',
            'Clear cell Meningioma',
            'Pleuropulmonary Blastoma',
            'Round Cell Sarcoma, NOS',
            'Myxoid Chondrosarcoma',
            'Dendritic Cell Sarcoma',
            'Lung Neuroendocrine Tumor',
            'Tubular Stomach Adenocarcinoma',
            'Testicular Mesothelioma',
            'Low-Grade Endometrial Stromal Sarcoma',
            'Uterine Mucinous Carcinoma',
            'Mixed Type Metaplastic Breast Cancer',
            'Villoglandular Adenocarcinoma of the Cervix',
            'Perivascular Epithelioid Cell Tumor',
            'Cancer of Unknown Primary, NOS',
            'Intimal Sarcoma',
            'Proximal-Type Epithelioid Sarcoma',
            'Pancreatoblastoma',
            'Signet Ring Mucinous Carcinoma',
            'Mixed Germ Cell Tumor',
            'Ampullary Carcinoma',
            'Breast Invasive Carcinosarcoma, NOS',
            'Congenital Nevus',
            'Adrenocortical Carcinoma',
            'Rosette-forming Glioneuronal Tumor of the Fourth Ventricle',
            'Glassy Cell Carcinoma of the Cervix',
            'Uterine Smooth Muscle Tumor of Uncertain Malignant Potential',
            'Large Cell Neuroendocrine Carcinoma',
            'Salivary Carcinoma, Other',
            'Hepatocellular Carcinoma',
            'High-Grade Neuroendocrine Carcinoma of the Ovary',
            'Neurofibroma',
            'Telangiectatic Osteosarcoma',
            'Well-Differentiated Liposarcoma',
            'Gliosarcoma',
            'Pilocytic Astrocytoma',
            'Salivary Duct Carcinoma',
            'Undifferentiated Uterine Sarcoma',
            'Pleura',
            'Dysgerminoma',
            'Spindle Cell Carcinoma of the Lung',
            'Peritoneal Mesothelioma',
            'Low-Grade Serous Ovarian Cancer',
            'Mucinous Borderline Ovarian Tumor',
            'Myxoma',
            'Endometrioid Carcinoma',
            'Lung Adenosquamous Carcinoma',
            'Gestational Trophoblastic Disease',
            'Malignant Lymphoma',
            'Clear Cell Ependymoma',
            'Polyembryoma',
            'Cervical Clear Cell Carcinoma',
            'Hepatoblastoma',
            'Small Cell Bladder Cancer',
            'Poorly Differentiated Vaginal Carcinoma',
            'Sarcomatoid Renal Cell Carcinoma',
            'Neuroendocrine Carcinoma, NOS',
            'Genitourinary Mucosal Melanoma',
            'Secondary Osteosarcoma',
            'Penis',
            'Ganglioglioma',
            'Prostate',
            'Cervical Adenoid Basal Carcinoma',
            'Adenoid Cystic Carcinoma of the Lung',
            'Aggressive Digital Papillary Adenocarcinoma',
            'Hepatocellular Carcinoma plus Intrahepatic Cholangiocarcinoma',
            'Papillary Thyroid Cancer',
            'Salivary Gland-Type Tumor of the Lung',
            'Colonic Type Adenocarcinoma of the Appendix',
            'Esophageal Squamous Cell Carcinoma',
            'Medulloblastoma with Extensive Nodularity',
            'Intestinal Ampullary Carcinoma',
            'Metaplastic Adenocarcinoma with Spindle Cell Differentiation',
            'Lung Carcinoid',
            'Chondroblastoma',
            'Germinoma',
            'Choroid Plexus Papilloma',
            'Ampulla of Vater',
            'Mucinous Ovarian Cancer',
            'Mucinous Adenocarcinoma of the Vulva/Vagina',
            'Larynx Squamous Cell Carcinoma',
            'High-Grade Glioma, NOS',
            'Urachal Carcinoma',
            'Interdigitating Dendritic Cell Sarcoma',
            'Uterine Myxoid Leiomyosarcoma',
            'Uveal Melanoma',
            'Esophagogastric Adenocarcinoma',
            'Ganglioneuroblastoma',
            'Cervical Adenoid Cystic Carcinoma',
            'Skin Adnexal Carcinoma',
            'Epithelioid Hemangioendothelioma',
            'Epithelial Type Metaplastic Breast Cancer',
            'Pleural Mesothelioma, Sarcomatoid Type',
            'Serous Ovarian Cancer',
            'Granulosa Cell Tumor',
            'Breast Carcinoma with Signet Ring',
            'Inflammatory Myofibroblastic Tumor',
            'Adenosquamous Carcinoma of the Pancreas',
            'Immature Teratoma',
            'Medulloblastoma',
            'Renal Small Cell Carcinoma',
            'Non-Seminomatous Germ Cell Tumor',
            'Papillary Tumor of the Pineal Region',
            'Adenoid Cystic Carcinoma',
            'Chordoid Meningioma',
            'Papillary Meningioma',
            'Pancreatic Adenocarcinoma',
            'Malignant Fibrothecoma',
            'Ovarian Carcinosarcoma/Malignant Mixed Mesodermal Tumor',
            'Signet Ring Cell Type of the Appendix',
            'Angiocentric Glioma',
            'Miscellaneous Neuroepithelial Tumor',
            'Thyroid',
            'Giant Cell Carcinoma of the Lung',
            'Large Cell/Anaplastic Medulloblastoma',
            'Metaplastic Adenosquamous Carcinoma',
            'Embryonal Tumor with Abundant Neuropil and True Rosettes',
            'Ewing Sarcoma',
            'Hypopharynx Squamous Cell Carcinoma',
            'Glomangiosarcoma',
            'Renal Non-Clear Cell Carcinoma',
            'High-Grade Neuroendocrine Carcinoma of the Colon and Rectum',
            'Inflammatory Myofibroblastic Lung Tumor',
            'Basaloid Penile Squamous Cell Carcinoma',
            'Soft Tissue',
            'Other Uterine Tumor',
            'Adrenocortical Adenoma',
            'Unclassified Renal Cell Carcinoma',
            'Brenner Tumor',
            'Mucoepidermoid Carcinoma of the Lung',
            'Wilms Tumor',
            'Uterine Epithelioid Leiomyosarcoma',
            'Pancreatic Neuroendocrine Tumor',
            'Peritoneum',
            'Desmoplastic Infantile Ganglioglioma',
            'Salivary Adenocarcinoma',
            'Perihilar Cholangiocarcinoma',
            'Ovarian Seromucinous Adenoma',
            'Poorly Differentiated Thyroid Cancer',
            'Pleomorphic Liposarcoma',
            'Endometrial Stromal Sarcoma',
            'Rhabdoid Meningioma',
            'Dedifferentiated Chondrosarcoma',
            'Pulmonary Lymphangiomyomatosis',
            'Small Cell Carcinoma of the Stomach',
            'Meningioma',
            'Breast Ductal Carcinoma In Situ',
            'Steroid Cell Tumor, NOS',
            'Adenoid Cystic Breast Cancer',
            'Serous Cystadenoma of the Pancreas',
            'Diffuse Glioma',
            'Merkel Cell Carcinoma',
            'Serous Borderline Ovarian Tumor',
            'Gastrointestinal Neuroendocrine Tumors',
            'Borderline Phyllodes Tumor of the Breast',
            'Thymic Neuroendocrine Tumor',
            'Sweat Gland Adenocarcinoma',
            'Basal Cell Carcinoma',
            'Neuroblastoma',
            'Aggressive Angiomyxoma',
            'Adenocarcinoma of the Gastroesophageal Junction',
            'Pleomorphic Rhabdomyosarcoma',
            'Renal Oncocytoma',
            'Lung Adenocarcinoma',
            'Sarcoma, NOS',
            'Clear Cell Ovarian Cancer',
            'Mesenchymal Chondrosarcoma',
            'Porocarcinoma/Spiroadenocarcinoma',
            'Skin',
            'Mucosal Melanoma of the Urethra',
            'Penile Squamous Cell Carcinoma',
            'Oligoastrocytoma',
            'Atypical Choroid Plexus Papilloma',
            'Anal Squamous Cell Carcinoma',
            'Pilomyxoid Astrocytoma',
            'Pancreatobiliary Ampullary Carcinoma',
            'High-Grade Neuroepithelial Tumor',
            'Urethral Cancer',
            'Anaplastic Ganglioglioma',
            'Poorly Differentiated Non-Small Cell Lung Cancer',
            'Germ Cell Tumor, Brain',
            'Peripheral Nervous System',
            'Extraskeletal Myxoid Chondrosarcoma',
            'Poroma/Acrospiroma',
            'Cervical Serous Carcinoma',
            'Kidney',
            'Colon Adenocarcinoma',
            'Mixed Ampullary Carcinoma',
            'Pleural Mesothelioma, Biphasic Type',
            'Epithelial-Myoepithelial Carcinoma',
            'Anaplastic Thyroid Cancer',
            'Retinoblastoma',
            'Adenosquamous Carcinoma of the Stomach',
            'Myopericytoma',
            'Intrahepatic Cholangiocarcinoma',
            'Mucosal Melanoma of the Vulva/Vagina',
            'Sinonasal Adenocarcinoma',
            'Bowel',
            'Mucoepidermoid Carcinoma',
            'Pituitary Adenoma',
            'Prostate Small Cell Carcinoma',
            'Alveolar Rhabdomyosarcoma',
            'Bladder Squamous Cell Carcinoma',
            'Mucinous Adenocarcinoma of the Appendix',
            'Pineal Parenchymal Tumor of Intermediate Differentiation',
            'Metaplastic Breast Cancer',
            'Brenner Tumor, Borderline',
            'Dysplastic Gangliocytoma of the Cerebellum/Lhermitte-Duclos Disease',
            'Pleural Mesothelioma',
            'Renal Clear Cell Carcinoma with Sarcomatoid Features',
            'Proliferating Pilar Cystic Tumor',
            'Biliary Tract',
            'Osteosarcoma',
            'Craniopharyngioma, Papillary Type',
            'Cervical Neuroendocrine Tumor',
            'Follicular Thyroid Cancer',
            'Appendiceal Adenocarcinoma',
            'High-Grade Surface Osteosarcoma',
            'Brenner Tumor, Malignant',
            'Rhabdomyosarcoma',
            'Mixed Cancer Types',
            'Olfactory Neuroblastoma',
            'Uterine Sarcoma, Other',
            'Tenosynovial Giant Cell Tumor Diffuse Type',
            'Mixed Ovarian Carcinoma',
            'Cervical Rhabdomyosarcoma',
            'Hurthle Cell Thyroid Cancer',
            'Diffuse Type Stomach Adenocarcinoma',
            'Uterine Neuroendocrine Carcinoma',
            'Head and Neck Mucosal Melanoma',
            'Medullomyoblastoma',
            'Undifferentiated Stomach Adenocarcinoma',
            'Hepatocellular Adenoma',
            'Hemangioma',
            'Uterine Mesonephric Carcinoma',
            'Intestinal Type Mucinous Carcinoma',
            'Ovarian Cancer, Other',
            'Sinonasal Undifferentiated Carcinoma',
            'Thymic Carcinoma',
            'Myoepithelial Carcinoma',
            'Uterine Clear Cell Carcinoma',
            'Glioblastoma',
            'Chondrosarcoma',
            'Small Intestinal Carcinoma',
            'Anaplastic Meningioma',
            'Gastrointestinal Stromal Tumor'
        ])
        self.trial_match_df = None

        # test match tree
        self.simple_mutation_match_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Mutation'
                    }
                },
                {
                    'clinical': {
                        s.mt_age: '>=18',
                        s.mt_diagnosis: 'Leiomyosarcoma'
                    }
                }
            ]
        }
        self.all_solid_tumor_match_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Mutation'
                    }
                },
                {
                    'clinical': {
                        s.mt_age: '>=18',
                        s.mt_diagnosis: '_SOLID_'
                    }
                }
            ]
        }
        self.all_male_match_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Mutation'
                    }
                },
                {
                    'clinical': {
                        s.mt_gender: 'Male',
                        s.mt_diagnosis: '_SOLID_'
                    }
                }
            ]
        }
        self.all_female_match_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Mutation'
                    }
                },
                {
                    'clinical': {
                        s.mt_gender: 'Female',
                        s.mt_diagnosis: '_LIQUID_'
                    }
                }
            ]
        }
        self.pediatric_cnv_match_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Copy Number Variation'
                    }
                },
                {
                    'clinical': {
                        s.mt_age: '<18',
                        s.mt_diagnosis: 'Hodgkin Lymphoma'
                    }
                }
            ]
        }
        self.exclusion_tree = {
            'and': [
                {
                    'or': [
                        {
                            'and': [
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'BRAF',
                                        s.mt_protein_change: 'p.V600E',
                                        s.mt_variant_category: '!Mutation'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'TP53',
                                        s.mt_protein_change: 'p.R278W'
                                    }
                                }
                            ]
                        },
                        {
                            'and': [
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'BRAF',
                                        s.mt_protein_change: 'p.V600D',
                                        s.mt_variant_category: '!Mutation'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'TP53',
                                        s.mt_protein_change: 'R278X'
                                    }
                                }
                            ]
                        }
                    ]
                },
                {
                    'and': [
                        {
                            'clinical': {
                                s.mt_age: '<18',
                                s.mt_diagnosis: '_SOLID_'
                            }
                        },
                        {
                            'clinical': {
                                s.mt_diagnosis: '!Melanoma'
                            }
                        }
                    ]
                }
            ]
        }

        self.complex_match_tree = {
            'and': [
                {
                    'and': [
                        {
                            'or': [
                                {
                                    'genomic': {
                                        'exon': 19,
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_variant_class: 'In_Frame_Del'
                                    }
                                },
                                {
                                    'genomic': {
                                        'exon': 19,
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_variant_class: 'In_Frame_Ins'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.E709K'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.G719S'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.G719A'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.G719D'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.G719C'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.L747P'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.A750P'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.S768I'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.L833V'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.L858R'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.L861R'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.L861R'
                                    }
                                }
                            ]
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'ALK',
                                s.mt_variant_category: '!Structural Variation'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'ROS1',
                                s.mt_variant_category: '!Structural Variation'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'RET',
                                s.mt_variant_category: '!Structural Variation'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'MET',
                                s.mt_cnv_call: 'High Amplification',
                                s.mt_variant_category: '!Copy Number Variation'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'MET',
                                s.mt_variant_category: '!Structural Variation'
                            }
                        },
                        {
                            'genomic': {
                                'exon': 14,
                                s.mt_hugo_symbol: 'MET',
                                s.mt_variant_category: '!Mutation',
                                s.mt_variant_class: 'In_Frame_Del'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'ERBB2',
                                s.mt_variant_category: '!Mutation',
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'BRAF',
                                s.mt_protein_change: 'p.V600E',
                                s.mt_variant_category: '!Mutation',
                            }
                        }
                    ]
                },
                {
                    'clinical': {
                        s.mt_age: '>=18',
                        s.mt_diagnosis: 'Lung Adenocarcinoma'
                    }
                }
            ]
        }

        self.simple_and_tree = {
            'and': [
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'EGFR',
                        s.mt_variant_category: 'Mutation',
                        s.mt_protein_change: 'p.E709K'
                    }
                },
                {
                    'genomic': {
                        s.mt_hugo_symbol: 'BRAF',
                        s.mt_variant_category: 'Mutation',
                        s.mt_protein_change: 'p.V600E'
                    }
                }
            ]
        }
        self.simple_or_tree = {
            'and': [
                {
                    'or': [
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'EGFR',
                                s.mt_variant_category: 'Mutation',
                                s.mt_protein_change: 'p.E709K'
                            }
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'BRAF',
                                s.mt_variant_category: 'Mutation',
                                s.mt_protein_change: 'p.V600E'
                            }
                        }
                    ]
                },
                {
                    'clinical': {
                        s.mt_diagnosis: 'Lung Adenocarcinoma'
                    }
                }
            ]

        }
        self.complex_or_tree = {
            'and': [
                {
                    'or': [
                        {
                            'and': [
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'EGFR',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.E709K'
                                    }
                                },
                                {
                                    'genomic': {
                                        s.mt_hugo_symbol: 'AKT',
                                        s.mt_variant_category: 'Mutation',
                                        s.mt_protein_change: 'p.E17K'
                                    }
                                }
                            ]
                        },
                        {
                            'genomic': {
                                s.mt_hugo_symbol: 'BRAF',
                                s.mt_variant_category: 'Mutation',
                                s.mt_protein_change: 'p.V600E'
                            }
                        }
                    ]
                },
                {
                    'clinical': {
                        s.mt_diagnosis: 'Lung Adenocarcinoma'
                    }
                }
            ]

        }

    def add_test_trials(self):
        """
        Adds all test trials to the database.

        :return: {null}
        """
        trial_docs = []
        for r, d, f in os.walk(self.data_path):
            for fn in f:
                if fn.endswith('yml'):
                    with open(os.path.join(r, fn), 'r') as ff:
                        trial_docs.append(yaml.load(ff))

        self.db.trial.insert_many(trial_docs)

    def add_test_trial_matches(self):
        """
        Add all test trials found in the test data path to this class

        :return: {null}
        """

        trial_match_data = {
            kn.tm_sample_id_col: ['TEST-SAMPLE-TM-1'],
            kn.tm_trial_protocol_no_col: ['00-000'],
            kn.tm_mrn_col: ['100100'],
            kn.tm_vital_status_col: [vital_status_allowed_vals[0]],
            kn.tm_trial_accrual_status_col: [s.match_accrual_status_open_val],
            kn.tm_sort_order_col: [0],
        }
        cols, data = trial_match_data.keys(), trial_match_data.values()
        self.trial_match_df = pd.DataFrame(data=[data], columns=cols)

    def load_trial(self, trial):
        """
        Load the trial at the given yml_path into this class

        :param yml_path: {trial} Trial name (e.g. 10-113)
        :return: {dict}
        """
        with open(os.path.join(self.data_path, '%s.yml' % trial), 'r') as f:
            doc = yaml.load(f)

        return doc

    def _find(self, query, proj=None, table='testSamples'):
        if proj is None:
            proj = self.proj

        return self.db[table].find_one(query, proj)

    def _findall(self, query, proj=None, table='testSamples'):
        if proj is None:
            proj = self.proj

        return list(self.db[table].find(query, proj))

    def _findalls(self, query, proj=None, table='testSamples'):
        return sorted([i[kn.sample_id_col] for i in self._findall(query=query, proj=proj, table=table)])

    @staticmethod
    def _print(query):
        print json.dumps(query, indent=4, default=str)


def get_demo_trial_matches():

    tm1 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '111-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation (variant-level match)
    tm2 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '222-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.E695Q',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 4,
                kn.allele_fraction_col: 0.38
            }
        ],
        kn.cnv_list_col: [
            {
                kn.hugo_symbol_col: 'BRAF',
                kn.variant_category_col: s.variant_category_cnv_val,
                kn.actionability_col: 'actionable'
            }
        ]
    }  # Tier 4 mutation with actionable CNV
    tm3 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '333-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.E695Q',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 4,
                kn.allele_fraction_col: 0.38
            }
        ],
        kn.cnv_list_col: [
            {
                kn.hugo_symbol_col: 'BRAF',
                kn.variant_category_col: s.variant_category_cnv_val,
            }
        ]
    }  # Tier 4 mutation with null-tier CNV
    tm4 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '444-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.G469A',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 2,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 2 mutation (variant-level match)
    tm5 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '555-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.N581S',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 3,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 3 mutation (variant-level match)
    tm6 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '666-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.E695Q',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 4,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 4 mutation (variant-level match)
    tm7 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '777-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'gene',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation (gene-level match)
    tm8 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '888-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: '_solid_',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation (all solid tumor match)
    tm9 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '999-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Massachussetts General Hospital',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation (MGH coordinating center)
    tm10 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '000-111',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation
    tm11 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '0001-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.mmr_status_col: s.mmr_status_deficient_val,
        kn.mutation_list_col: [
            {
                kn.mr_reason_level_col: 'variant',
                kn.hugo_symbol_col: 'BRAF',
                kn.protein_change_col: 'p.V600E',
                kn.variant_category_col: s.variant_category_mutation_val,
                kn.variant_class_col: 'Missense_Mutation',
                kn.transcript_exon_col: 15,
                kn.chromosome_col: '1',
                kn.strand_col: '-',
                kn.position_col: 140453136,
                kn.ref_allele_col: 'A',
                kn.cdna_change_col: 'c.1799T>A',
                kn.tier_col: 1,
                kn.allele_fraction_col: 0.38
            }
        ]
    }  # Tier 1 mutation with MMR Deficient signature
    tm12 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '0002-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.wt_genes_col: [
            {
                kn.hugo_symbol_col: 'BRAF',
            }
        ]
    }  # BRAF wild-type
    tm13 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '0003-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female',
        kn.mr_coordinating_center_col: 'Dana-Farber Cancer Institute',
        kn.sv_list_col: [
            {
                kn.sv_comment_col: 'BRAF',
            }
        ]
    }  # BRAF structural variation
    tm14 = {
        kn.sample_id_col: '111',
        kn.mrn_col: '111',
        kn.tm_trial_protocol_no_col: '0004-000',
        kn.mr_trial_level_col: 'step',
        kn.mr_trial_step_code_col: '1',
        kn.tm_trial_accrual_status_col: 'open',
        kn.vital_status_col: 'alive',
        kn.oncotree_primary_diagnosis_name_col: 'Cutaneous Melanoma',
        kn.mr_diagnosis_level_col: 'specific',
        kn.gender_col: 'Female'
    }  # clinical-only match
    trial_matches = [tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8, tm9, tm10, tm11, tm12, tm13, tm14]
    return trial_matches
