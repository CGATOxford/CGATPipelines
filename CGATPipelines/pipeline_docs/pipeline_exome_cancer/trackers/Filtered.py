import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class Snp(ExomeTracker):

    pattern = "(\S*)_mutect_snp_annotated$"

    def __call__(self, track, slice=None):

        tables = self.getValues("SELECT tbl_name FROM sqlite_master")
        tables = list(set([x for x in tables if re.match(".*_annotations", x)]))

        if len(tables) > 0:
            sql_columns = ["%s.%s" % (x, x.replace("_annotations", ""))
                           for x in tables]
            annotations_select_cmd = "%s," % ",".join(sql_columns)

            sql_joins = ["%s ON A.SNPEFF_GENE_NAME = %s.gene_id" % (x, x)
                         for x in tables]
            annotations_join_cmd = "LEFT JOIN %s" % " LEFT JOIN ".join(sql_joins)

        else:
            annotations_select_cmd = ""
            annotations_join_cmd = ""

        statement = '''
        SELECT A.CHROM AS Chr, A.POS AS Pos,
        A.SNPEFF_GENE_NAME AS Gene,
        A.SNPEFF_EXON_ID AS Exon,
        A.REF, A.ALT,
        A.SNPEFF_IMPACT AS Impact, A.SNPEFF_GENE_BIOTYPE AS Biotype,
        SNPEFF_AMINO_ACID_CHANGE AS AA_change,
        SNPEFF_CODON_CHANGE AS Codon_change,
        %(annotations_select_cmd)s
        C.type as NCG, C.cancer_type, D.*,
        B.n_ref_count AS Normal_Ref, B.n_alt_count AS Normal_Alt,
        B.t_ref_count AS Tumor_Ref, B.t_alt_count AS Tumor_Alt
        FROM %(track)s_mutect_snp_annotated AS A
        JOIN %(track)s_call_stats AS B
        ON A.CHROM = B.contig AND A.POS = B.position
        LEFT OUTER JOIN cancergenes as C
        ON A.SNPEFF_GENE_NAME = C.symbol
        LEFT OUTER JOIN eBio_studies_gene_frequencies as D
        ON A.SNPEFF_GENE_NAME = D.gene
        %(annotations_join_cmd)s
        ''' % locals()

        return self.getAll(statement)


class Indel(ExomeTracker):

    pattern = "(\S*)_indels_annotated$"

    def __call__(self, track, slice=None):

        tables = self.getValues("SELECT tbl_name FROM sqlite_master")
        tables = list(set([x for x in tables if re.match(".*_annotations", x)]))

        if len(tables) > 0:
            sql_columns = ["%s.%s" % (x, x.replace("_annotations", ""))
                           for x in tables]
            annotations_select_cmd = "%s," % ",".join(sql_columns)

            sql_joins = ["%s ON A.SNPEFF_GENE_NAME = %s.gene_id" % (x, x)
                         for x in tables]
            annotations_join_cmd = "LEFT JOIN %s" % " LEFT JOIN ".join(sql_joins)

        else:
            annotations_select_cmd = ""
            annotations_join_cmd = ""

        statement = '''
        SELECT A.CHROM AS Chr, A.POS AS Pos,
        A.SNPEFF_GENE_NAME AS Gene,
        A.SNPEFF_EXON_ID AS Exon,
        A.REF, A.ALT,
        A.SNPEFF_IMPACT AS Impact, A.SNPEFF_GENE_BIOTYPE AS Biotype,
        A.SNPEFF_AMINO_ACID_CHANGE AS AA_change,
        A.SNPEFF_CODON_CHANGE AS Codon_change,
        %(annotations_select_cmd)s
        B.type as NCG, B.cancer_type,  C.*,
        A.NORMAL_DP AS Normal_depth,
        A.TUMOR_DP AS Tumor_depth,
        A.NORMAL_TAR as Normal_Ref, A.NORMAL_TIR as Normal_Alt,
        A.TUMOR_TAR as Tumor_Ref, A.TUMOR_TIR as Tumor_Alt
        FROM %(track)s_indels_annotated AS A
        LEFT OUTER JOIN cancergenes as B
        ON A.SNPEFF_GENE_NAME = B.symbol
        LEFT OUTER JOIN eBio_studies_gene_frequencies as C
        ON A.SNPEFF_GENE_NAME = C.gene
        %(annotations_join_cmd)s
        ''' % locals()

        return self.getAll(statement)


class FilterSummary(ExomeTracker):

    pattern = "(\S*)_mutect_filtering_summary$"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT justification, IFNULL(single, 0.0) AS Single,
        IFNULL(combination, 0.0) AS Combination
        FROM %(track)s_mutect_filtering_summary
        ;
        ''' % locals()

        return self.getAll(statement)
