from src.utilities import settings as s
from src.data_store import key_names as kn

# mutation data
mutation_missense_data = {
    kn.variant_category_col: s.variant_category_mutation_val,
    kn.hugo_symbol_col: 'BRAF',
    kn.chromosome_col: '7',
    kn.position_col: 140453135,
    kn.strand_col: '-',
    kn.transcript_exon_col: 15,
    kn.variant_class_col: s.variant_class_missense_mutation_val,
    kn.protein_change_col: 'p.V600E',
    kn.cdna_change_col: 'c.1799_1800TG>AA',
    kn.cdna_transcript_id_col: 'NM_004333',
    kn.alt_allele_col: 'TT',
    kn.ref_allele_col: 'CA',
    kn.allele_fraction_col: 0.11,
    kn.transcript_src_col: 'UCSC',
    kn.coverage_col: 100,
    kn.somatic_status_col: 'Likely Somatic',
    kn.tier_col: 1,
    kn.entrez_id_col: '673',
    kn.is_best_effect_col: False,
    'extraCol': None
}

mutation_nonsense_data = {
    kn.variant_category_col: s.variant_category_mutation_val,
    kn.hugo_symbol_col: 'BRAF',
    kn.chromosome_col: '7',
    kn.position_col: 140482826,
    kn.strand_col: '-',
    kn.transcript_exon_col: 9,
    kn.variant_class_col: s.variant_class_nonsense_mutation_val,
    kn.protein_change_col: 'p.R437*',
    kn.cdna_change_col: 'c.1309C>T',
    kn.cdna_transcript_id_col: 'NM_004333',
    kn.alt_allele_col: 'A',
    kn.ref_allele_col: 'G',
    kn.allele_fraction_col: 0.11,
    kn.transcript_src_col: 'UCSC',
    kn.coverage_col: 100,
    kn.somatic_status_col: 'Likely Somatic',
    kn.tier_col: 1,
    kn.entrez_id_col: '673',
    kn.is_best_effect_col: False,
    'extraCol': None
}

# cnv data
cnv_heterozygous_del_data = {
    kn.variant_category_col: s.variant_category_cnv_val,
    kn.hugo_symbol_col: 'WHSC1',
    kn.cytoband_col: '4 p',
    kn.cnv_call_col: 'Heterozygous deletion',
    kn.cnv_band_col: '4p16.3',
    kn.copy_count_col: None,
    kn.cnv_row_id_col: 5130,
    kn.actionability_col: 'investigational',
    'extraCol': None
}

# sv data
sv_data = {
    kn.variant_category_col: s.variant_category_sv_val,
    kn.sv_comment_col: 'NTRK1 fusion detected'
}

# signature data
signature_mmr_deficient_data = {
    kn.mmr_status_col: 'Deficient (MMR-D / MSI-H)',
    kn.tobacco_status_col: 'Yes',
    kn.tmz_status_col: None,
    kn.pole_status_col: 'TEST'
}

signature_mmr_none_data = {
    kn.mmr_status_col: None
}

# wildtype data
wt1_data = {kn.hugo_symbol_col: 'BRAF'}
wt2_data = {kn.hugo_symbol_col: 'EGFR'}

# low coverage data
pertinent_negative_data = {
    kn.hugo_symbol_col: 'BRAF',
    kn.coverage_type_col: 'PN',
    kn.roi_type_col: 'C',
    kn.panel_col: 'GLIOMA',
    kn.codon_col: 466,
    kn.show_codon_col: True,
    kn.lc_coverage_col: None,
    kn.low_coverage_exon_col: 11,
}
pertinent_negative_v2_data = {
    kn.hugo_symbol_col: 'EGFR',
    kn.coverage_type_col: 'PN',
    kn.roi_type_col: 'C',
    kn.panel_col: 'THORACIC',
    kn.codon_col: 709,
    kn.show_codon_col: True,
    kn.lc_coverage_col: None,
    kn.low_coverage_exon_col: 18,
}

pertinent_undercovered_data = {
    kn.hugo_symbol_col: 'PTEN',
    kn.coverage_type_col: 'PLC',
    kn.roi_type_col: 'G',
    kn.panel_col: 'GLIOMA',
    kn.codon_col: None,
    kn.show_exon_col: True,
    kn.lc_coverage_col: 11.6923,
    kn.low_coverage_exon_col: 3,
}

additional_undercovered_data = {
    kn.hugo_symbol_col: 'PTK2B',
    kn.coverage_type_col: 'NPLC',
    kn.roi_type_col: 'G',
    kn.panel_col: 'GLIOMA',
    kn.codon_col: None,
    kn.show_exon_col: True,
    kn.lc_coverage_col: 11.6923,
    kn.low_coverage_exon_col: 18,
}
