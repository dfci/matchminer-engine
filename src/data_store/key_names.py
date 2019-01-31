# ------------------ #
# Genomic key names  #
# ------------------ #

# Mutations key names
variant_category_col = 'variantCategory'
hugo_symbol_col = 'hugoSymbol'
chromosome_col = 'chromosome'
position_col = 'position'
strand_col = 'strand'
transcript_exon_col = 'transcriptExon'
variant_class_col = 'variantClassification'
protein_change_col = 'proteinChange'
cdna_change_col = 'cDNAChange'
cdna_transcript_id_col = 'cDNATranscriptId'
alt_allele_col = 'alternateAllele'
ref_allele_col = 'referenceAllele'
ref_residue_col = 'referenceResidue'
allele_fraction_col = 'alleleFraction'
transcript_src_col = 'transcriptSource'
coverage_col = 'coverage'
somatic_status_col = 'somaticStatus'
tier_col = 'tier'
entrez_id_col = 'entrezId'
is_best_effect_col = 'isBestEffect'

# CNVs key names
cytoband_col = 'cytoband'
cnv_call_col = 'cnvCall'
cnv_band_col = 'cnvBand'
copy_count_col = 'copyCount'
cnv_row_id_col = 'cnvRowId'
actionability_col = 'actionability'

# SVs key names
sv_comment_col = 'structuralVariantComment'
sv_type_col = 'svType'
sv_lp_gene = 'lpHugoSymbol'
sv_lp_chromosome_col = 'lpChromosome'
sv_lp_strand = 'lpStrand'
sv_lp_intergenic = 'lpIntergenic'
sv_lp_position = 'lpPosition'
sv_lp_band = 'lpBand'
sv_lp_type = 'lpType'
sv_lp_type_number = 'lpTypeNumber'
sv_rp_gene = 'rpHugoSymbol'
sv_rp_chromosome_col = 'rpChromosome'
sv_rp_strand = 'rpStrand'
sv_rp_intergenic = 'rpIntergenic'
sv_rp_position = 'rpPosition'
sv_rp_band = 'rpBand'
sv_rp_type = 'rpType'
sv_rp_type_number = 'rpTypeNumber'


# Low Coverage key names
coverage_type_col = 'coverageType'
roi_type_col = 'roiType'
panel_col = 'panel'
codon_col = 'codon'
show_codon_col = 'showCodon'
low_coverage_exon_col = 'exon'
show_exon_col = 'showExon'
lc_coverage_col = 'coverageLC'

# ------------------ #
# Clinical key names #
# ------------------ #

# patient identification information
mrn_col = 'mrn'
sample_id_col = 'sampleId'
block_no_col = 'blockNumber'
data_push_id_col = 'dataPushId'
alt_mrn_col = 'altMRN'
patient_id_col = 'patientId'
powerpath_patient_id_col = 'powerpathPatientId'

# patient personal information
gender_col = 'gender'
first_name_col = 'firstName'
last_name_col = 'lastName'
birth_date_col = 'birthDate'
vital_status_col = 'vitalStatus'
first_last_col = 'firstLast'
last_first_col = 'lastFirst'

# patient appointment information
disease_center_descr_col = 'diseaseCenterDescr'
ord_physician_name_col = 'ordPhysicianName'
ord_physician_npi_col = 'ordPhysicianNPI'
ord_physicians_email_col = 'ordPhysicianEmail'
pathologist_name_col = 'pathologistName'
last_visit_date_col = 'lastVisitDate'

# patient oncopanel information
panel_version_col = 'panelVersion'
report_comment_col = 'reportComment'
date_received_at_seq_center_col = 'dateReceivedAtSeqCenter'
report_date_col = 'reportDate'
report_version_col = 'reportVersion'
test_type_col = 'testType'
qc_result_col = 'qcResult'
cnv_results_col = 'cnvResults'
snv_results_col = 'snvResults'
cancer_type_percentile_col = 'cancerTypePercentile'
all_profile_percentile_col = 'allProfilePercentile'
metamain_count_col = 'metamainCount'
case_count_col = 'caseCount'
pdf_layout_version_col = 'pdfLayoutVersion'

# sequencing information
ref_genomic_col = 'referenceGenome'
total_aligned_reads_col = 'totalAlignedReads'
pct_target_base_col = 'pctTargetBase'
mean_sample_coverage_col = 'meanSampleCoverage'
total_reads_col = 'totalReads'
tmb_col = 'tumorMutationalBurdenPerMegabase'

# patient clinical information
oncotree_primary_diagnosis_col = 'oncotreePrimaryDiagnosis'
oncotree_primary_diagnosis_name_col = 'oncotreePrimaryDiagnosisName'
tumor_purity_percent_col = 'tumorPurityPercent'
oncotree_biopsy_site_col = 'oncotreeBiopsySite'
oncotree_biopsy_site_type_col = 'oncotreeBiopsySiteType'
oncotree_biopsy_site_name_col = 'oncotreeBiopsySiteName'
oncotree_biopsy_site_meta_col = 'oncotreeBiopsySiteMeta'
oncotree_biopsy_site_color_col = 'oncotreeBiopsySiteColor'
oncotree_primary_diagnosis_meta_col = 'oncotreePrimaryDiagnosisMeta'
oncotree_primary_diagnosis_color_col = 'oncotreePrimaryDiagnosisColor'

# patient consent information
q1_consent_col = 'question1Consent'
q2_consent_col = 'question2Consent'
q3_consent_col = 'question3Consent'
q4_consent_col = 'question4Consent'
q5_consent_col = 'question5Consent'
cris_consent_col = 'crisConsent'
consent_17000_col = 'consent17000'

# ------------------ #
# Samples key names  #
# ------------------ #
# wild-type information
wt_genes_col = 'wildTypeGenes'

# variant information
mutation_list_col = 'mutations'
cnv_list_col = 'cnvs'
sv_list_col = 'svs'
pertinent_negatives_list_col = 'pertinent_negatives'
pertinent_undercovered_list_col = 'pertinent_undercovered'
additional_undercovered_list_col = 'additional_undercovered'

# mutational signature information
mmr_status_col = 'mmrStatus'
ms_status_col = 'msStatus'
tobacco_status_col = 'tobaccoStatus'
tmz_status_col = 'tmzStatus'
pole_status_col = 'polEStatus'
apobec_status_col = 'apobecStatus'
uva_status_col = 'uvaStatus'
all_signatures = [
    mmr_status_col, ms_status_col, tobacco_status_col, tmz_status_col, pole_status_col, apobec_status_col,
    uva_status_col
]

# ----------------------- #
# Trial matches key names #
# ----------------------- #
tm_sample_id_col = 'sampleId'
tm_trial_protocol_no_col = 'trialProtocolNo'
tm_mrn_col = 'mrn'
tm_vital_status_col = 'vitalStatus'
tm_trial_accrual_status_col = 'trialAccrualStatus'
tm_sort_order_col = 'sortOrder'
tm_match_reasons_col = 'matchReasons'

# ----------------------- #
# Match reasons key names #
# ----------------------- #
mr_trial_level_col = 'trialLevel'
mr_trial_step_code_col = 'trialStepCode'
mr_trial_arm_code_col = 'trialArmCode'
mr_trial_dose_code_col = 'trialDoseCode'
mr_low_coverage_list_col = 'lowCoverage'
mr_diagnosis_col = 'diagnosis'
mr_age_col = 'age'
mr_gender_col = 'gender'
mr_reason_level_col = 'level'
mr_diagnosis_level_col = 'diagnosisLevel'
mr_coordinating_center_col = 'coordinatingCenter'
genomic_exclusion_reasons_col = 'genomic_exclusion_reasons'
clinical_exclusion_reasons_col = 'clinical_exclusion_reasons'
