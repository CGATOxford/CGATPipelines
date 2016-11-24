import pandas as pd
import numpy as np
import sqlite3

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *

###############################################################################
# parse params
###############################################################################
DATABASE = P.get('', P.get('sql_backend', 'sqlite:///./csvdb'))
ANNOTATIONS_DATABASE = P.get('annotations_database')
###############################################################################
# trackers
###############################################################################


class SleuthResults(IsoformTracker):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL"

    def __call__(self, track, slice=None):

        quantifier = track.split("_")[-2]

        statement = '''
        SELECT A.gene_name, A.gene_id, A.transcript_id, B.reason AS flagged,
        A.control_mean AS expression, A.fold, A.l2fold, A.p_value,
        A.p_value_adj, A.significant, A.transcript_biotype
        FROM %(track)s_DEresults AS A
        LEFT JOIN %(quantifier)s_flagged_transcripts AS B
        ON A.transcript_id = B.transcript_id
        %(where)s
        ORDER BY A.significant DESC, A.l2fold ASC
        '''

        return self.getAll(statement)


class SleuthResultsSig(SleuthResults):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL AND significant == 1 AND ABS(l2fold) > 1"


class SleuthAll(IsoformTracker):

    pattern = ""
    table = ""

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.*, B.reason as flagged
        FROM all_%(table)s_%(track)s AS A
        LEFT JOIN %(track)s_flagged_transcripts AS B
        ON A.transcript_id = B.transcript_id
        '''
        return self.getAll(statement)


class SleuthCountsAll(SleuthAll):
    pattern = "all_counts_(.*)"
    table = "counts"


class SleuthTpmAll(SleuthAll):
    pattern = "all_tpm_(.*)"
    table = "tpm"


class SummarisedResults(IsoformTracker):

    pattern = "(.*)_DEresults$"

    def __call__(self, track, slice=None):

        quantifier = track.split("_")[-2]
        design = "_".join(track.split("_")[:-2]) + "_design"

        select_results = '''SELECT A.*,
        B.p_value, B.p_value_adj, B.l2fold, B.transcript_biotype
        FROM all_tpm_%(quantifier)s AS A
        LEFT JOIN %(track)s_DEresults AS B
        ON A.transcript_id = B.transcript_id
        WHERE B.p_value<0.05'''

        results_df = pd.DataFrame(self.getAll(select_results))

        select_design = '''SELECT * FROM %(design)s'''
        design_df = self.getDataFrame(select_design)

        for group in set(design_df['_group']):
            group_tracks = design_df[design_df["_group"] == group]['track']
            group_tracks = [x.replace("-", "_") for x in group_tracks]

            results_df["group_%s_mean" % group] = results_df[group_tracks].mean(axis=1)
            results_df["group_%s_stdev" % group] = results_df[group_tracks].std(axis=1)

        return results_df


class SleuthAllGenes(IsoformTracker):

    pattern = ""
    table = ""

    def __call__(self, track, slice=None):

        statement = "SELECT * FROM all_gene_expression_%(table)s_%(track)s"

        return self.getDataFrame(statement)


class SleuthCountsAllGenes(SleuthAllGenes):
    pattern = "all_gene_expression_counts_(.*)"
    table = "counts"


class SleuthTpmAllGenes(SleuthAllGenes):
    pattern = "all_gene_expression_tpm_(.*)"
    table = "tpm"


class Deseq2Results(IsoformTracker):

    pattern = "(.*)_deseq2_DE_results$"
    direction = ""
    where = "WHERE p_value NOT NULL"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT
        gene_name, gene_id,
        control_mean AS expression,
        fold, l2fold,
        p_value, p_value_adj,
        significant
        FROM %(track)s_deseq2_DE_results
        WHERE p_value NOT NULL
        ORDER BY significant DESC, l2fold ASC
        '''

        return self.getAll(statement)


class Deseq2ResultsSig(Deseq2Results):

    pattern = "(.*)_deseq2_DE_results$"
    direction = ""
    where = "WHERE p_value NOT NULL AND significant == 1 AND ABS(l2fold) > 1"
