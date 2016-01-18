import pandas as pd
import numpy as np

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *

###############################################################################
# parse params
###############################################################################
DATABASE = P.get('', P.get('sql_backend', 'sqlite:///./csvdb'))

###############################################################################
# trackers
###############################################################################

class SleuthResults(IsoformTracker):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.gene_name, A.gene_id, A.transcript_id, B.reason AS flagged,
        A.control_mean AS expression, A.fold, A.l2fold, A.p_value,
        A.p_value_adj, A.significant, A.transcript_biotype
        FROM %(track)s_DEresults AS A
        LEFT JOIN kallisto_flagged_transcripts AS B
        ON A.test_id = B.transcript_id
        %(where)s
        ORDER BY A.significant DESC, A.l2fold ASC
        '''

        return self.getAll(statement)


class SleuthResultsSig(SleuthResults):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL AND significant == 1 AND ABS(l2fold) > 1"


class SleuthAll(IsoformTracker):

    table = ""

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.*, B.reason as flagged
        FROM %(table)s AS A
        LEFT JOIN kallisto_flagged_transcripts AS B
        ON A.transcript_id = B.transcript_id
        '''
        return self.getAll(statement)


class SleuthCountsAll(SleuthAll):
    table = "all_counts"


class SleuthTpmAll(SleuthAll):
    table = "all_tpm"


class SummarisedResults(IsoformTracker):

    pattern = "(.*)_DEresults$"

    def __call__(self, track, slice=None):

        select_results = '''SELECT A.*,
        B.p_value, B.p_value_adj, B.l2fold, B.transcript_biotype
        FROM %(track)s_tpm AS A
        LEFT JOIN %(track)s_DEresults AS B
        ON A.transcript_id = B.transcript_id
        WHERE B.p_value<0.05'''

        results_df = pd.DataFrame(self.getAll(select_results))

        select_design = '''SELECT * FROM %(track)s_design'''
        design_df = self.getDataFrame(select_design)

        for group in set(design_df['_group']):
            group_tracks = design_df[design_df["_group"] == group]['track']
            group_tracks = [x.replace("-", "_") for x in group_tracks]

            results_df["group_%s_mean" % group] = results_df[group_tracks].mean(axis=1)
            results_df["group_%s_stdev" % group] = results_df[group_tracks].std(axis=1)

        return results_df


#class GeneLevelExpression(IsoformTracker):

#    pattern = "(.*)_design$"

#    def __call__(self, track, slice=None):

#        select_results = '''SELECT * FROM %(track)s_tpm''' % locals()

#        results_df = pd.DataFrame(self.getAll(select_results))

#        grouped_df = results_df.groupby(["gene_id", "gene_name"])
#        grouped_df = pd.DataFrame(grouped_df.aggregate(sum))

#        select_design = '''SELECT * FROM %(track)s_design'''
#        design_df = self.getDataFrame(select_design)

#        for group in set(design_df['_group']):
#            group_tracks = design_df[design_df["_group"]==group]['track']
#            group_tracks = [x.replace("-", "_") for x in group_tracks]

#            grouped_df["group_%s_mean" % group] = grouped_df[group_tracks].mean(axis=1)
#            grouped_df["group_%s_stdev" % group] = grouped_df[group_tracks].std(axis=1)

#        return grouped_df
