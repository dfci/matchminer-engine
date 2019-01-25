import datetime as dt

from src.utilities import settings as s
from src.data_store import key_names as kn
from tests.test_load_service import mutation_missense_data, mutation_nonsense_data, cnv_heterozygous_del_data, \
    sv_data, pertinent_undercovered_data, pertinent_negative_data, pertinent_negative_v2_data, \
    additional_undercovered_data


valid_mutation_missense_data = mutation_missense_data.copy()
valid_mutation_nonsense_data = mutation_nonsense_data.copy()
valid_mutation_missense_data[kn.ref_residue_col] = 'p.V600'
valid_mutation_missense_data[kn.ref_residue_col] = None
del valid_mutation_missense_data['extraCol']
del valid_mutation_nonsense_data['extraCol']

valid_cnv_data = cnv_heterozygous_del_data.copy()
del valid_cnv_data['extraCol']

samples_data = {

    # patient identification information
    kn.mrn_col: '123',
    kn.sample_id_col: 'DEV-SAMPLE-01',
    kn.block_no_col: 'DEV-BLOCK-01',
    kn.data_push_id_col: '2000-01-01 12:00:00',
    kn.alt_mrn_col: None,
    kn.patient_id_col: None,
    kn.powerpath_patient_id_col: None,

    # patient personal information
    kn.gender_col: 'Male',
    kn.first_name_col: 'FirstNameDev',
    kn.last_name_col: 'LastNameDev',
    kn.birth_date_col: dt.datetime(year=2000, day=1, month=1, hour=4),
    kn.vital_status_col: 'alive',

    # patient appointment information
    kn.disease_center_descr_col: 'GYNECOLOGIC ONCOLOGY',
    kn.ord_physician_name_col: 'OrdLastNameDev OrdFirstNameDev A. M.D.',
    kn.ord_physician_npi_col: 1000000000,
    kn.ord_physicians_email_col: 'ordemaildev@partners.org',
    kn.pathologist_name_col: 'PathFirstNameDev A. PathLastNameDev M.D.',
    kn.last_visit_date_col: None,

    # patient oncopanel information
    kn.panel_version_col: 3,
    kn.report_comment_col: 'DEV Report Comment',
    kn.date_received_at_seq_center_col: dt.datetime(year=2000, day=1, month=1, hour=5),
    kn.report_date_col: dt.datetime(year=2000, day=1, month=1, hour=5),
    kn.report_version_col: 2,
    kn.test_type_col: 'PROFILECOHORT',
    kn.qc_result_col: 'PASS',
    kn.cnv_results_col: 'Has CNV Results',
    kn.snv_results_col: 'Has SNV Results',
    kn.cancer_type_percentile_col: 0.1,
    kn.all_profile_percentile_col: 0.1,
    kn.metamain_count_col: 100,
    kn.case_count_col: 1000,
    kn.pdf_layout_version_col: 2,

    # sequencing information
    kn.ref_genomic_col: 'hg19',
    kn.total_aligned_reads_col: 1000000,
    kn.pct_target_base_col: 0.1,
    kn.mean_sample_coverage_col: 100,
    kn.total_reads_col: 10000000,
    kn.tmb_col: 1,

    # patient clinical information
    kn.oncotree_primary_diagnosis_col: 'LMS',
    kn.oncotree_primary_diagnosis_name_col: 'Leiomyosarcoma',
    kn.tumor_purity_percent_col: 0.1,
    kn.oncotree_biopsy_site_col: 'UTERUS',
    kn.oncotree_biopsy_site_type_col: 'Unspecified',
    kn.oncotree_biopsy_site_name_col: 'Uterus',
    kn.oncotree_biopsy_site_meta_col: None,
    kn.oncotree_biopsy_site_color_col: None,
    kn.oncotree_primary_diagnosis_meta_col: 'Soft Tissue Sarcoma',
    kn.oncotree_primary_diagnosis_color_col: None,

    # patient consent information
    kn.q1_consent_col: 'Y',
    kn.q2_consent_col: 'Y',
    kn.q3_consent_col: 'Y',
    kn.q4_consent_col: 'Y',
    kn.q5_consent_col: 'Y',
    kn.cris_consent_col: 'Y',
    kn.consent_17000_col: None,

    # wild-type information
    kn.wt_genes_col: [{
        kn.variant_category_col: s.variant_category_wt_val,
        kn.hugo_symbol_col: 'BRAF',
    }],

    # variant information
    kn.mutation_list_col: [valid_mutation_nonsense_data, valid_mutation_missense_data],
    kn.cnv_list_col: [valid_cnv_data],
    kn.sv_list_col: [sv_data],
    kn.pertinent_negatives_list_col: [pertinent_negative_data, pertinent_negative_v2_data],
    kn.pertinent_undercovered_list_col: [pertinent_undercovered_data],
    kn.additional_undercovered_list_col: [additional_undercovered_data],

    # mutational signature information
    kn.mmr_status_col: 'Deficient',
    kn.ms_status_col: 'MSI-H',
    kn.tobacco_status_col: 'Yes',
    kn.tmz_status_col: None,
    kn.pole_status_col: None,
    kn.apobec_status_col: None,
    kn.uva_status_col: None
}
