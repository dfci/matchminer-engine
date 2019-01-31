from src.utilities import settings as s
from src.data_store import key_names as kn
from src.data_store.samples_data_model import misc_signature_allowed_vals


def assess_inclusion(node_value):
    """
    Assess if the given node value sholud be treated as an inclusion or an exclusion criterion.
    If it begins with "!" treat as an exclusion criterion. Otherwise, treat as an inclusion criterion.

    :param node_value: {str}
    :return: {bool} True if inclusion, False if exclusion
    """
    return not node_value.startswith('!')


def sanitize_exclusion_vals(val):
    """
    Remove the preceding "!" from exclusion criteria values.

    :param val: {str}
    :return: {str}
    """
    if val.startswith('!'):
        return val[1:]
    else:
        return val


def normalize_variant_category_val(val):
    """
    Normalize the variant category value to what is expected in the samples table in the database.

    :param val: {str}
    :return: {str}
    """
    variant_category_dict = {
        s.mt_cnv_val: s.variant_category_cnv_val,
        s.mt_mut_val: s.variant_category_mutation_val,
        s.mt_sv_val: s.variant_category_sv_val,
        s.mt_any_vc_val: s.variant_category_any_val
    }
    return variant_category_dict[val]


def normalize_cnv_call_val(val):
    """
    Normalize the cnv call value to what is expected in the samples table in the database.

    :param val: {str}
    :return: {str}
    """
    cnv_call_dict = {
        s.mt_high_amp_val: s.cnv_call_high_amp,
        s.mt_homo_del_val: s.cnv_call_homo_del,
        s.mt_hetero_del_val: s.cnv_call_hetero_del,
        s.mt_low_amp_val: s.cnv_call_gain
    }
    return cnv_call_dict[val]


def normalize_signature_vals(signature_type, signature_val):
    """
    Normalize the mutational signature type and value to what is expected in the samples table in the database.

    :param signature_type: {str} (e.g. mmr_status, ms_status, etc.)
    :param signature_val: {str} (e.g. (MMR-Deficient, MSI-H, etc.)
    :return: {tuple of str} (type, val)
    """
    signature_type_dict = {
        s.mt_mmr_status: kn.mmr_status_col,
        s.mt_ms_status: kn.ms_status_col,
        s.mt_tobacco_status: kn.tobacco_status_col,
        s.mt_tmz_status: kn.tmz_status_col,
        s.mt_pole_status: kn.pole_status_col,
        s.mt_apobec_status: kn.apobec_status_col,
        s.mt_uva_status: kn.uva_status_col
    }
    signature_val_dict = {
        s.mt_mmr_deficient_val: s.mmr_status_deficient_val,
        s.mt_mmr_proficient_val: s.mmr_status_proficient_val,
        s.mt_msi_high_val: s.ms_status_msih_val,
        s.mt_mss_val: s.ms_status_mss_val,
    }
    for val in misc_signature_allowed_vals:
        signature_val_dict[val] = val

    return signature_type_dict[signature_type], signature_val_dict[signature_val]


def add_unwind(node, include):
    """
    Adds the unwind command for the aggregation pipeline for the given digraph node.

    :param node: {Directed graph node}
    :param include: {bool}
    :return: {null}
    """
    if 'unwind' not in node:
        node['unwind'] = '$%s' % node['match_reason'].keys()[0].split('.')[0] \
            if 'match_reason' in node and node['match_reason'] is not None else None
    node['include'] = include
    return node


def reformat_matches(matches):
    """Reformat matches so that they fit cerberus validation."""
    variant_keys = [kn.mutation_list_col, kn.cnv_list_col, kn.sv_list_col, kn.wt_genes_col]
    for match in matches:
        for variant_key in variant_keys:
            if variant_key in match:
                match[variant_key] = [match[variant_key]]

    return matches
