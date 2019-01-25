from src.utilities import settings as s
from src.data_store import key_names as kn

# shared rules
variant_category_rules = {
    'type': 'string',
    'required': True,
    'allowed': [
        s.variant_category_mutation_val,
        s.variant_category_cnv_val,
        s.variant_category_sv_val,
        s.variant_category_wt_val
    ]
}
vital_status_allowed_vals = ['alive', 'deceased']
misc_signature_allowed_vals = ['Yes', 'No', 'Cannot assess', 'insufficient variants', None]
tier_allowed_vals = [1, 2, 3, 4]
actionability_allowed_vals = ['actionable', 'investigational', None]
low_coverage_allowed_vals = [s.pertinent_negative_val, s.pertinent_low_coverage_val, s.additional_low_coverage_val]
gender_allowed_vals = ['Male', 'Female']

# Genomic schemas
mutations_schema = {
    kn.variant_category_col: variant_category_rules,
    kn.hugo_symbol_col: {'type': 'string', 'required': True},
    kn.chromosome_col: {'type': 'string', 'required': False},
    kn.position_col: {'type': 'integer', 'required': False},
    kn.strand_col: {'type': 'string', 'required': True, 'allowed': ['+', '-', None], 'nullable': True},
    kn.transcript_exon_col: {'type': 'integer', 'required': True},
    kn.variant_class_col: {'type': 'string', 'required': True},
    kn.protein_change_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.cdna_change_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.cdna_transcript_id_col: {'type': 'string', 'nullable': True},
    kn.alt_allele_col: {'type': 'string', 'required': True},
    kn.ref_allele_col: {'type': 'string', 'required': True},
    kn.ref_residue_col: {'type': 'string', 'nullable': True},
    kn.allele_fraction_col: {'type': 'float', 'required': True, 'nullable': True},
    kn.transcript_src_col: {'type': 'string', 'nullable': True},
    kn.coverage_col: {'type': 'integer', 'nullable': True},
    kn.somatic_status_col: {'type': 'string', 'nullable': True},
    kn.tier_col: {'type': 'integer', 'allowed': tier_allowed_vals, 'required': True},
    kn.entrez_id_col: {'type': 'string', 'nullable': True},
    kn.is_best_effect_col: {'type': 'boolean', 'required': True}
}

cnvs_schema = {
    kn.variant_category_col: variant_category_rules,
    kn.hugo_symbol_col: {'type': 'string', 'required': True},
    kn.cytoband_col: {'type': 'string', 'required': False, 'nullable': True},
    kn.cnv_call_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.cnv_band_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.copy_count_col: {'type': 'integer', 'required': False, 'nullable': True},
    kn.cnv_row_id_col: {'type': 'integer', 'required': False, 'nullable': True},
    kn.actionability_col: {'type': 'string', 'required': False, 'nullable': True, 'allowed': actionability_allowed_vals}
}

svs_schema = {
    kn.variant_category_col: variant_category_rules,
    kn.sv_comment_col: {'type': 'string', 'required': False},
    kn.sv_type_col: {'type': 'string', 'allowed': ['Deletion', 'rearrangement', 'indel']},
    kn.sv_lp_gene: {'type': 'string', 'required': False},
    kn.sv_lp_chromosome_col: {'type': 'string'},
    kn.sv_lp_strand: {'type': 'string', 'allowed': ['+', '-', None], 'nullable': True},
    kn.sv_lp_intergenic: {'type': 'boolean', 'nullable': True},
    kn.sv_lp_position: {'type': 'integer'},
    kn.sv_lp_band: {'type': 'string'},
    kn.sv_lp_type: {'type': 'string', 'allowed': ['intron', 'exon', None], 'nullable': True},
    kn.sv_lp_type_number: {'type': 'integer'},
    kn.sv_rp_gene: {'type': 'string', 'required': False},
    kn.sv_rp_chromosome_col: {'type': 'string'},
    kn.sv_rp_strand: {'type': 'string', 'allowed': ['+', '-', None], 'nullable': True},
    kn.sv_rp_intergenic: {'type': 'boolean', 'nullable': True},
    kn.sv_rp_position: {'type': 'integer'},
    kn.sv_rp_band: {'type': 'string'},
    kn.sv_rp_type: {'type': 'string', 'allowed': ['intron', 'exon', None], 'nullable': True},
    kn.sv_rp_type_number: {'type': 'integer'},

}
wts_schema = {
    kn.variant_category_col: variant_category_rules,
    kn.hugo_symbol_col: {'type': 'string', 'required': True},
}

low_coverage_schema = {
    kn.hugo_symbol_col: {'type': 'string', 'required': True},
    kn.coverage_type_col: {'type': 'string', 'required': True, 'allowed': low_coverage_allowed_vals},
    kn.roi_type_col: {'type': 'string', 'allowed': ['C', 'G', 'R', 'E', 'M', None], 'nullable': True, 'required': True},
    kn.panel_col: {'type': 'string', 'required': True},
    kn.codon_col: {'type': 'integer', 'nullable': True},
    kn.show_codon_col: {'type': 'boolean', 'nullable': True},
    kn.show_exon_col: {'type': 'boolean', 'nullable': True},
    kn.lc_coverage_col: {'type': 'float', 'nullable': True},
    kn.low_coverage_exon_col: {'type': 'integer', 'nullable': True},
}

# Samples schema
samples_schema = {

    # patient identification information
    kn.mrn_col: {'type': 'string', 'required': True},
    kn.sample_id_col: {'type': 'string', 'required': True},
    kn.block_no_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.data_push_id_col: {'type': 'string', 'required': False, 'nullable': True},
    kn.alt_mrn_col: {'type': 'string', 'nullable': True},
    kn.patient_id_col: {'type': 'string', 'nullable': True},
    kn.powerpath_patient_id_col: {'type': 'string', 'nullable': True},

    # patient personal information
    kn.gender_col: {'type': 'string', 'required': True, 'allowed': gender_allowed_vals},
    kn.first_name_col: {'type': 'string', 'required': True},
    kn.last_name_col: {'type': 'string', 'required': True},
    kn.birth_date_col: {'type': 'datetime', 'required': True},
    kn.vital_status_col: {'type': 'string', 'required': True, 'allowed': vital_status_allowed_vals},
    kn.first_last_col: {'type': 'string', 'readonly': True},
    kn.last_first_col: {'type': 'string', 'readonly': True},

    # patient appointment information
    kn.disease_center_descr_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.ord_physician_name_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.ord_physician_npi_col: {'type': 'integer', 'required': True, 'nullable': True},
    kn.ord_physicians_email_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.pathologist_name_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.last_visit_date_col: {'type': 'datetime', 'nullable': True},

    # patient oncopanel information
    kn.panel_version_col: {'type': 'integer', 'required': True},
    kn.report_comment_col: {'type': 'string', 'required': True, 'nullable': True},
    kn.date_received_at_seq_center_col: {'type': 'datetime', 'required': True},
    kn.report_date_col: {'type': 'datetime', 'required': True},
    kn.report_version_col: {'type': 'integer', 'required': True, 'nullable': True},
    kn.test_type_col: {'type': 'string', 'required': True},
    kn.qc_result_col: {'type': 'string', 'required': True},
    kn.cnv_results_col: {'type': 'string'},
    kn.snv_results_col: {'type': 'string'},
    kn.cancer_type_percentile_col: {'type': 'float', 'nullable': True},
    kn.all_profile_percentile_col: {'type': 'float', 'nullable': True},
    kn.metamain_count_col: {'type': 'integer', 'nullable': True},
    kn.case_count_col: {'type': 'integer', 'nullable': True},
    kn.pdf_layout_version_col: {'type': 'integer', 'allowed': [1, 2]},

    # sequencing information
    kn.ref_genomic_col: {'type': 'string'},
    kn.total_aligned_reads_col: {'type': 'integer', 'required': True},
    kn.pct_target_base_col: {'type': 'float', 'required': True},
    kn.mean_sample_coverage_col: {'type': 'integer', 'required': True},
    kn.total_reads_col: {'type': 'integer', 'required': True},
    kn.tmb_col: {'type': 'float', 'nullable': True},

    # patient clinical information
    kn.oncotree_primary_diagnosis_col: {'type': 'string', 'required': True},
    kn.oncotree_primary_diagnosis_name_col: {'type': 'string', 'required': True},
    kn.tumor_purity_percent_col: {'type': 'float', 'required': True},
    kn.oncotree_biopsy_site_col: {'type': 'string', 'nullable': True},
    kn.oncotree_biopsy_site_type_col: {'type': 'string', 'nullable': True},
    kn.oncotree_biopsy_site_name_col: {'type': 'string', 'nullable': True},
    kn.oncotree_biopsy_site_meta_col: {'type': 'string', 'nullable': True},
    kn.oncotree_biopsy_site_color_col: {'type': 'string', 'nullable': True},
    kn.oncotree_primary_diagnosis_meta_col: {'type': 'string', 'nullable': True},
    kn.oncotree_primary_diagnosis_color_col: {'type': 'string', 'nullable': True},

    # patient consent information
    kn.q1_consent_col: {'type': 'string', 'required': True},
    kn.q2_consent_col: {'type': 'string', 'required': True},
    kn.q3_consent_col: {'type': 'string', 'required': True},
    kn.q4_consent_col: {'type': 'string', 'required': True},
    kn.q5_consent_col: {'type': 'string', 'required': False},
    kn.cris_consent_col: {'type': 'string', 'required': True},
    kn.consent_17000_col: {'type': 'string', 'nullable': True},

    # wild-type information
    kn.wt_genes_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': wts_schema}},

    # variant information
    kn.mutation_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': mutations_schema}},
    kn.cnv_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': cnvs_schema}},
    kn.sv_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': svs_schema}},
    kn.pertinent_negatives_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': low_coverage_schema}},
    kn.pertinent_undercovered_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': low_coverage_schema}},
    kn.additional_undercovered_list_col: {'type': 'list', 'schema': {'type': 'dict', 'schema': low_coverage_schema}},

    # mutational signature information
    kn.mmr_status_col: {
        'type': 'string',
        'allowed': [
            s.mmr_status_cannot_assess_val,
            s.mmr_status_indeterminate_val,
            s.mmr_status_proficient_val,
            s.mmr_status_deficient_val,
            None
        ],
        'nullable': True
    },
    kn.ms_status_col: {
        'type': 'string',
        'allowed': [
            s.mmr_status_cannot_assess_val,
            s.mmr_status_indeterminate_val,
            s.ms_status_mss_val,
            s.ms_status_msih_val,
            None
        ],
        'nullable': True
    },
    kn.tobacco_status_col: {'type': 'string', 'allowed': misc_signature_allowed_vals, 'nullable': True},
    kn.tmz_status_col: {'type': 'string', 'allowed': misc_signature_allowed_vals, 'nullable': True},
    kn.pole_status_col: {'type': 'string', 'allowed': misc_signature_allowed_vals, 'nullable': True},
    kn.apobec_status_col: {'type': 'string', 'allowed': misc_signature_allowed_vals, 'nullable': True},
    kn.uva_status_col: {'type': 'string', 'allowed': misc_signature_allowed_vals, 'nullable': True}
}
