from src.utilities import settings as s
from src.data_store import key_names as kn
from src.services.match_service.query_utils.clinical_utils import ClinicalUtils
from src.services.match_service.query_utils.genomic_utils import GenomicUtils


class ProjUtils(ClinicalUtils, GenomicUtils):

    def __init__(self):
        ClinicalUtils.__init__(self)
        GenomicUtils.__init__(self)

        self.clinical_inclusion_dict = {
            True: self.create_clinical_inclusion_proj,
            False: self.create_clinical_exclusion_proj
        }
        self.clinical_proj_dict = {
            s.mt_diagnosis: self.diagnosis_key,
            s.mt_age: self.age_key,
            s.mt_gender: self.gender_key
        }

        self.genomic_inclusion_dict = {
            True: self.create_genomic_inclusion_proj,
            False: self.create_genomic_exclusion_proj
        }
        self.genomic_proj_dict = {
            s.mt_hugo_symbol: self.hugo_symbol_key
        }

        self.proj = {
            '_id': 0,
            kn.sample_id_col: 1,
            kn.mrn_col: 1,
            kn.vital_status_col: 1
        }
        self.exclusion_key = 'ne'

    def create_clinical_proj(self, include=True, **kwargs):
        """
        Return either inclusive or exclusive clinical projetion.

        :param include {bool}
        :param kwargs
        :return: {dict}
        """
        subproj = self.clinical_inclusion_dict[include](**kwargs)
        if not include:
            return subproj

        proj = self.proj.copy()
        for k, v in subproj.iteritems():
            proj[k] = v

        return proj

    def create_clinical_inclusion_proj(self, **kwargs):
        """
        Create MongoDB projection to return only the matching diagnosis from the sample record.

        :param kwargs:
            - keys {list of str}

        :return: {dict}
        """
        return {self.clinical_proj_dict[k]: 1 for k in kwargs['keys']}

    def create_clinical_exclusion_proj(self, **kwargs):
        """
        Return indicator that matching records matched because they lacked the specified criteria

        :param kwargs:
            - keys {list of str}
            - vals {list of any type}

        :return: {dict}
        """
        return {self.clinical_proj_dict[k]: v for k, v in zip(kwargs['keys'], kwargs['vals'])}

    def create_genomic_proj(self, include=True, **kwargs):
        """
        Return either inclusive or exclusive genomic projetion.

        :param include: {bool}
        :param kwargs:
        :return: {dict}
        """
        subproj, match_reasons = self.genomic_inclusion_dict[include](**kwargs)
        if not include:
            return subproj, match_reasons

        proj = self.proj.copy()
        for k, v in subproj.iteritems():
            proj[k] = v

        return proj, match_reasons

    @staticmethod
    def create_genomic_inclusion_proj(**kwargs):
        """
        Create MongoDB projection to return only the matching genomic variants from the sample record.

        :param kwargs:
            - query {dict}

        :return: {dict}
        """
        variant_category = kwargs['query'].keys()[0]
        match_reason = {
            '%s.%s' % (variant_category, k): v for k, v in kwargs['query'][variant_category]['$elemMatch'].iteritems()
        }

        proj_keys = [kn.hugo_symbol_col, kn.protein_change_col, kn.variant_class_col, kn.cnv_call_col, kn.ref_residue_col]
        proj = {
            '%s.%s' % (variant_category, k): 1 for k in proj_keys
        }
        return proj, match_reason

    def create_genomic_exclusion_proj(self, **kwargs):
        """
        Return indicator that matching records matched because they lacked the specified criteria

        :param kwargs:
            - keys {list of str}
            - vals {list of any type}

        :return: {dict}
        """

        proj = {}
        for k, v in zip(kwargs['keys'], kwargs['vals']):
            if v is None:
                continue
            elif k in self.genomic_proj_dict:
                proj[self.genomic_proj_dict[k]] = v
            else:
                if k in [kn.mutation_list_col, kn.cnv_list_col, kn.sv_list_col, kn.wt_genes_col]:
                    k = kn.variant_category_col

                proj[k] = v

        return proj, None

    def create_any_variant_proj(self, gene_name, include=True, query=None):
        """
        Create an inclusive or exclusive genomic projection for an any variant category query

        :param gene_name: {str}
        :param include: {bool}
        :param query: {dict or null}
        :return: {dict}
        """
        any_variant_dict = {True: self.create_genomic_proj, False: self.create_any_variant_exclusion_proj}
        return any_variant_dict[include](include=include, gene_name=gene_name, query=query)

    def create_any_variant_exclusion_proj(self, gene_name, **kwargs):
        """
        Create an exclusive genomic projection for an exclusion query against any variant category.

        :param gene_name: {str}
        :return: {dict}
        """
        mutation_proj = self.create_genomic_exclusion_proj(keys=[kn.mutation_list_col, kn.hugo_symbol_col],
                                                           vals=[s.variant_category_mutation_val, gene_name])
        cnv_proj = self.create_genomic_exclusion_proj(keys=[kn.cnv_list_col, kn.hugo_symbol_col],
                                                      vals=[s.variant_category_cnv_val, gene_name])
        return {'$and': [mutation_proj, cnv_proj]}
