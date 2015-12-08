import pandas as pd
import numpy as np

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *


class imagesTracker(TrackerImages):

    '''Convience Tracker for globbing images for gallery plot'''
    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]


class TranscriptBiotypeSummary(IsoformTracker):

    pattern = "(.*)_DEresults$"
    direction = ""

    def __call__(self, track, slice=None):

        statement = '''SELECT transcript_id, transcript_biotype, significant
                       FROM %(track)s_DEresults
                       WHERE l2fold %(direction)s 0;'''

        df = pd.DataFrame(self.getAll(statement))

        keep_biotypes = ["protein_coding", "retained_intron",
                         "processed_transcript", "nonsense_mediated_decay"]

        df = df[[x in keep_biotypes for x in df['transcript_biotype']]]

        grouped = df.groupby(['significant', 'transcript_biotype'])
        df_agg = grouped.aggregate({"transcript_id": 'count'})
        df_agg.columns = ["Count"]

        # TS: must be able to do this more succinctly!
        fraction = []
        for l in df_agg.groupby(level=0).agg(
                lambda x: list(x/float(np.sum(x)))).values.flatten():
            fraction.extend(l)

        cumsum_values = []
        for l in df_agg.groupby(level=0).agg(
                lambda x: list(np.cumsum(x)/float(np.sum(x)))).values.flatten():
            cumsum_values.extend(l)

        previous = 0
        cumsum_centred = []
        for value in cumsum_values:
            if previous == 1:
                previous = 0
            cumsum_centred.append((previous+value)/2)
            previous = value

        df_agg['fraction'] = fraction
        df_agg['cumsum'] = cumsum_values
        df_agg['cumsum_centres'] = cumsum_centred

        df_agg.reset_index(inplace=True)
        df_agg['significant'] = ["Significant" if x == 1
                                 else "Not significant"
                                 for x in df_agg['significant']]

        return df_agg


class TranscriptBiotypeSummaryUp(TranscriptBiotypeSummary):

    direction = ">"


class TranscriptBiotypeSummaryDown(TranscriptBiotypeSummary):

    direction = "<"


class TranscriptExpressionOrdered(IsoformTracker):

    pattern = "(.*)_tpm$"

    def __call__(self, track, slice=None):

        def ordered_log(array):
            array = sorted(array, reverse=True)
            return np.log10(array)

        statement = '''SELECT * FROM %(track)s_tpm;'''

        df = pd.DataFrame(self.getAll(statement))

        df.drop(["gene_id", "gene_name",
                 "transcript_biotype", "transcript_id"],
                axis=1, inplace=True)

        df = df.apply(ordered_log, axis=0)

        df["index"] = df.index

        df = pd.melt(df, id_vars=["index"])

        df = df.replace([np.inf, -np.inf], np.nan).dropna()
        df.index = ['all', ] * len(df)

        return df


class TranscriptNumberSamplesExpressed(IsoformTracker):

    pattern = "(.*)_tpm$"

    def __call__(self, track, slice=None):

        statement = '''SELECT * FROM %(track)s_tpm;'''

        df = pd.DataFrame(self.getAll(statement))

        df = df.set_index(["transcript_id"])
        df.drop(["gene_id", "gene_name", "transcript_biotype"],
                axis=1, inplace=True)

        final_df = pd.DataFrame()

        for threshold in (0.01, 0.1, 0.5, 1, 2, 10):
            df_tmp = pd.DataFrame(
                {"Count": df.apply(func=lambda row:
                                   sum([x > threshold for x in row]), axis=1),
                 "No_transcripts": range(0, len(df.index))})

            df_tmp = df_tmp.ix[df_tmp['Count'] > 0, :]

            df_tmp = df_tmp.groupby(["Count"]).count()
            df_tmp.reset_index(inplace=True)
            df_tmp["threshold"] = [threshold]*len(df_tmp.index)

            final_df = pd.concat([final_df, df_tmp])

        return final_df
