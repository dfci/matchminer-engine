import os
import sys
import json
import logging

from src.data_store import key_names as kn

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )

MONGO_URI = ""
MONGO_DBNAME = ""

# pull secrets
file_path = os.getenv("SECRETS_JSON", None)

if file_path is None:
    logging.error("ENVAR SECRETS_JSON not set")

if file_path is not None:
    with open(file_path) as fin:
        vars = json.load(fin)
        for name, value in vars.iteritems():
            setattr(sys.modules[__name__], name, value)


sample_collection_name = 'samples'
variant_category_mutation_val = 'MUTATION'
variant_category_cnv_val = 'CNV'
variant_category_sv_val = 'SV'
variant_class_missense_mutation_val = 'Missense_Mutation'
variant_class_nonsense_mutation_val = 'Nonsense_Mutation'
mmr_status_cannot_assess_val = 'Cannot assess'
mmr_status_indeterminate_val = 'Indeterminate (see note)'
mmr_status_proficient_val = 'Proficient'
mmr_status_deficient_val = 'Deficient'
ms_status_mss_val = 'MSS'
ms_status_msih_val = 'MSI-H'
pertinent_negative_val = 'PN'
pertinent_low_coverage_val = 'PLC'
additional_low_coverage_val = 'NPLC'
signature_cols = [
    kn.mmr_status_col,
    kn.ms_status_col,
    kn.tobacco_status_col,
    kn.tmz_status_col,
    kn.pole_status_col,
    kn.apobec_status_col,
    kn.uva_status_col
]