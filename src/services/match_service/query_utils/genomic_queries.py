import re

from src.utilities import settings as s
from src.data_store import key_names as kn
from src.services.match_service.query_utils.query_utils import QueryUtils
from src.services.match_service.query_utils.genomic_utils import GenomicUtils


class GenomicQueries(QueryUtils, GenomicUtils):

    def __init__(self):
        QueryUtils.__init__(self)
        GenomicUtils.__init__(self)

        self.inclusion_dict = {
            True: self.create_inclusion_query,
            False: self.create_exclusion_query
        }
        self.variant_level_query_dict = {
            True: self.create_variant_level_inclusion_query,
            False: self.create_variant_level_exclusion_query
        }
        self.exon_query_dict = {
            True: self.create_exon_inclusion_query,
            False: self.create_exon_exclusion_query
        }
        self.variant_class_query_dict = {
            True: self.create_variant_class_inclusion_query,
            False: self.create_variant_class_exclusion_query
        }

        self.variant_category_dict = {
            s.variant_category_mutation_val: kn.mutation_list_col,
            s.variant_category_wildcard_mutation_val: kn.mutation_list_col,
            s.variant_category_exon_val: kn.mutation_list_col,
            s.variant_category_variant_class_val: kn.mutation_list_col,
            s.variant_category_cnv_val: kn.cnv_list_col,
            s.variant_category_sv_val: kn.sv_list_col,
            s.variant_category_wt_val: kn.wt_genes_col
        }
        self.variant_type_col_dict = {
            s.variant_category_mutation_val: self.protein_change_key,
            s.variant_category_wildcard_mutation_val: self.ref_residue_key,
            s.variant_category_exon_val: self.transcript_exon_key,
            s.variant_category_cnv_val: self.cnv_call_key
        }

    def create_gene_level_query(self, gene_name, variant_category, include=True):
        """
        Create MongoDB query to find records by gene name and variant category

        :param gene_name: {str}
        :param variant_category: {str} (MUTATION, CNV, SV, WT)
        :param include: {bool}
        :return: {dict}
        """
        return self.inclusion_dict[include](variant_category=self.variant_category_dict[variant_category],
                                            key=self.hugo_symbol_key,
                                            val=gene_name)

    def create_variant_level_inclusion_query(self, variant_category, gene_name, variant_val):
        """
        Create MongoDB query that matches the specific variant specified.

        :param variant_category: {str} (e.g. MUTATION, EXON, WILDCARD_MUTATION, CNV)
        :param gene_name: {str}
        :param variant_val: {str}
        :return: {dict}
        """
        return {
            self.variant_category_dict[variant_category]: {
                '$elemMatch': {
                    self.hugo_symbol_key: gene_name,
                    self.variant_type_col_dict[variant_category]: variant_val
                }
            }
        }

    def create_variant_level_exclusion_query(self, variant_category, gene_name, variant_val):
        """
        Create MongoDB query that matches the specific variant specified.

        :param variant_category: {str} (e.g. MUTATION, WILDCARD_MUTATION, CNV)
        :param gene_name: {str}
        :param variant_val: {str}
        :return: {dict}
        """
        exclude_query = {
            self.variant_category_dict[variant_category]: {
                '$not': {
                    '$elemMatch': {
                        self.hugo_symbol_key: {'$eq': gene_name},
                        self.variant_type_col_dict[variant_category]: {'$eq': variant_val}
                    }
                }
            }
        }
        query = {'$or': [
            self.create_gene_level_query(gene_name=gene_name,
                                         variant_category=variant_category,
                                         include=False),
            exclude_query
        ]}

        return query

    def create_mutation_query(self, gene_name, protein_change, include=True):
        """
        Create MongoDB query to find SNV records by gene name and protein change

        :param gene_name: {str}
        :param protein_change: {str}
        :param include: {bool}
        :return: {dict}
        """
        return self.variant_level_query_dict[include](variant_category=s.variant_category_mutation_val,
                                                      gene_name=gene_name,
                                                      variant_val=protein_change)

    def create_wildcard_query(self, gene_name, protein_change, include=True):
        """
        Create MongoDB query to find samples that match the given wildcard protein change value

        :param gene_name: {str}
        :param protein_change: {str}
        :param include: {bool}
        :return: {dict}
        """
        return self.variant_level_query_dict[include](variant_category=s.variant_category_wildcard_mutation_val,
                                                      gene_name=gene_name,
                                                      variant_val=protein_change)

    def create_variant_class_query(self, gene_name, variant_class, include=True):
        """
        Create MongoDB query to find samples with mutations in the given gene and variant classification

        :param gene_name: {str}
        :param variant_class: {str}
        :param include: {bool}
        :return: {dict}
        """
        return self.variant_class_query_dict[include](gene_name=gene_name,
                                                      variant_class=variant_class)

    def create_variant_class_inclusion_query(self, gene_name, variant_class):
        """
        Create MongoDB query that matches the specific gene and variant classification.

        :param gene_name: {str}
        :param variant_class: {str}
        :return: {dict}
        """
        return {
            kn.mutation_list_col: {
                '$elemMatch': {
                    self.hugo_symbol_key: gene_name,
                    self.variant_class_key: variant_class
                }
            }
        }

    def create_variant_class_exclusion_query(self, gene_name, variant_class):
        """
        Create MongoDB query that matches the specific variant specified.

        :param gene_name: {str}
        :param exon: {int}
        :param variant_class: {str or null}
        :return: {dict}
        """
        exclude_query = {
            kn.mutation_list_col: {
                '$elemMatch': {
                    self.hugo_symbol_key: {'$eq': gene_name},
                    self.variant_class_key: {'$ne': variant_class}
                }
            }
        }
        query = {'$or': [
            self.create_gene_level_query(gene_name=gene_name,
                                         variant_category=s.variant_category_variant_class_val,
                                         include=False),
            exclude_query
        ]}
        return query

    def create_exon_query(self, gene_name, exon, variant_class=None, include=True):
        """
        Create MongoDB query to find samples with mutations in the given exon.

        :param gene_name: {str}
        :param exon: {int}
        :param variant_class: {str or null}
        :param include: {bool}
        :return: {dict}
        """
        return self.exon_query_dict[include](gene_name=gene_name,
                                             exon=exon,
                                             variant_class=variant_class)

    def create_exon_inclusion_query(self, gene_name, exon, variant_class=None):
        """
        Create MongoDB query that matches the specific variant specified.

        :param gene_name: {str}
        :param exon: {int}
        :param variant_class: {str or null}
        :return: {dict}
        """
        query = {
            kn.mutation_list_col: {
                '$elemMatch': {
                    self.hugo_symbol_key: gene_name,
                    self.transcript_exon_key: exon
                }
            }
        }
        if variant_class is not None:
            query[kn.mutation_list_col]['$elemMatch'][self.variant_class_key] = variant_class

        return query

    def create_exon_exclusion_query(self, gene_name, exon, variant_class=None):
        """
        Create MongoDB query that matches the specific variant specified.

        :param gene_name: {str}
        :param exon: {int}
        :param variant_class: {str or null}
        :return: {dict}
        """
        exclude_query = {
            kn.mutation_list_col: {
                '$elemMatch': {
                    self.hugo_symbol_key: {'$eq': gene_name},
                    self.transcript_exon_key: {'$ne': exon}
                }
            }
        }
        query = {'$or': [
            self.create_gene_level_query(gene_name=gene_name,
                                         variant_category=s.variant_category_exon_val,
                                         include=False),
            exclude_query
        ]}
        if variant_class is not None:
            query['$or'].append({
                kn.mutation_list_col: {
                    '$elemMatch': {
                        self.hugo_symbol_key: {'$eq': gene_name},
                        self.transcript_exon_key: {'$eq': exon},
                        self.variant_class_key: {'$ne': variant_class}
                    }
                }
            })

        return query

    def create_cnv_query(self, gene_name, cnv_call, include=True):
        """
        Create MongoDB query to find CNV records by gene name and, if supplied, cnv call

        :param gene_name: {str}
        :param cnv_call: {str}
        :param include: {bool}
        :return: {dict}
        """
        return self.variant_level_query_dict[include](variant_category=s.variant_category_cnv_val,
                                                      gene_name=gene_name,
                                                      variant_val=cnv_call)

    def create_sv_query(self, gene_name, include=True):
        """
        Create MongoDB query to find SV records by regex string matching of the given gene name
        to a free text structrual variant comment field.

        :param gene_name: {str}
        :param include: {bool}
        :return: {dict}
        """
        gene_name_regex = self.regex_compile_gene_name(gene_name)
        return self.inclusion_dict[include](variant_category=kn.sv_list_col,
                                            key=self.sv_comment_key,
                                            val=gene_name_regex)

    @staticmethod
    def create_mutational_signature_query(signature_type, signature_val):
        """
        Create MongoDB query to find records by the given mutational signature value.

        :param signature_type: {str} (e.g. MMR Status, TMZ Status, etc.)
        :param signature_val: {str}
        :return: {dict}
        """
        return {signature_type: signature_val}

    def create_low_coverage_query(self, gene_name, low_coverage_type, include=True):
        """
        Create MongoDB query to find samples with genes identified as low coverage.

        :param gene_name: {str}
        :param low_coverage_type: {str} (e.g. Pertinent Low Coverage, Additional Low Coverage)
        :param include: {bool}
        :return: {dict}
        """
        raise NotImplementedError

    def create_any_variant_query(self, gene_name, include=True):
        """
        Create MongoDB query to find samples with any variation in the specified gene.

        :param gene_name: {str}
        :param include: {bool}
        :return: {dict}
        """
        include_dict = {True: '$or', False: '$and'}
        return {include_dict[include]: [
            self.create_gene_level_query(gene_name=gene_name,
                                         variant_category=s.variant_category_mutation_val,
                                         include=include),
            self.create_gene_level_query(gene_name=gene_name,
                                         variant_category=s.variant_category_cnv_val,
                                         include=include),
        ]}
