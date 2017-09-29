import numpy as np
import pandas as pd
import random
import re
import itertools
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS
from CGATReport.ResultBlock import ResultBlocks, ResultBlock

DATABASE = PARAMS.get('scrnaseqqc_backend', PARAMS.get('sql_backend', 'sqlite:///./csvdb'))


class ScQCTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def getMetaData(self):
        '''
        Retrieve the experimental data table
        Create the sample names from the
        seqRun, Plate and Well columns
        '''

        meta_state = ''''
        SELECT * FROM experimental_design;
        '''

        meta_df = self.getDataFrame(meta_state)
        meta_df["track"] = ["_%s_%s_%s" % (meta_df.loc[x, "seqRun"],
                                           meta_df.loc[x, "Plate"],
                                           meta_df.loc[x, "Well"]) for x in meta_df.index]
        meta_df.index = meta_df["track"]

        return meta_df

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class ExpressionTracker(ScQCTracker):

    '''
    Tracker to pull out expression measurements
    '''

    def __init__(self, *args, **kwargs):
        # table_regex var is to be set for each sub class
        ScQCTracker.__init__(self, *args, **kwargs)


    def getExpressionTable(self, table_regex):
        '''
        retrieve the appropriate SQLite table
        '''

        table_state = """SELECT name FROM
        sqlite_master WHERE type='table';"""

        tables = self.getAll(table_state)

        tab_reg = re.compile(self.table_regex)
        # getAll returns an orderedDict with values as a single
        # tuple
        tracker_table = [tx for tx in tables.values()[0] if re.search(tab_reg,
                                                                      tx)][0]
        return tracker_table

    def __call__(self, track, slice=None):


        tracker_table = self.getExpressionTable(self.table_regex)
        statement = '''
        SELECT * FROM %(tracker_table)s;
        '''
        df = self.getDataFrame(statement)

        return df
        

class SailfishExpressionTable(ExpressionTracker):

    table_regex = r"_sailfish$"


class FeatureCountsExpressionTable(ExpressionTracker):

    table_regex = r"_feature_counts"


class SailfishSpikeInExpressionTable(SailfishExpressionTable):

    def getSpikeIns(self, dataframe):
        '''
        Extract just the spike in expression
        from a sailfish expression table
        '''

        spike_regex = re.compile("ERCC")
        spike_indx = [sx for sx in dataframe["Name"] if re.search(spike_regex,
                                                                  sx)]
        dataframe.index = dataframe["Name"]
        out_df = dataframe.loc[spike_indx]
        out_df.drop(labels="Name", inplace=True,
                    axis=1)

        return out_df

    def __call__(self, track, slice=None):

        tracker_table = self.getExpressionTable(self.table_regex)
        statement = '''
        SELECT * FROM %(tracker_table)s;
        '''

        df = self.getDataFrame(statement)

        fin_df = self.getSpikeIns(df)

        return fin_df

class FeatureCountsSpikeInExpressionTable(FeatureCountsExpressionTable):

    def getSpikeIns(self, dataframe):
        '''
        Extract just the spike in expression
        from a sailfish expression table
        '''

        spike_regex = re.compile("ERCC")
        spike_indx = [sx for sx in dataframe["gene_id"] if re.search(spike_regex,
                                                                     sx)]
        dataframe.index = dataframe["gene_id"]
        out_df = dataframe.loc[spike_indx]
        out_df.index = [ix.replace(".", "-") for ix in out_df.index]

        out_df.drop(labels="gene_id", inplace=True,
                    axis=1)

        return out_df

    def __call__(self, track, slice=None):

        tracker_table = self.getExpressionTable(self.table_regex)
        statement = '''
        SELECT * FROM %(tracker_table)s;
        '''

        df = self.getDataFrame(statement)

        fin_df = self.getSpikeIns(df)

        return fin_df


class CompareSpikeInExpression(ExpressionTracker):

    def __init__(self, *args, **kwargs):
        ScQCTracker.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):

        sailfish = SailfishSpikeInExpressionTable()
        feature = FeatureCountsSpikeInExpressionTable()

        sfish_df = sailfish(track)
        feat_df = feature(track)
        sfish_df["group"] = "sailfish"
        sfish_df["gene_id"] = sfish_df.index
        feat_df["group"] = "featureCounts"
        feat_df["gene_id"] = feat_df.index

        sailfish_melt = pd.melt(sfish_df,
                                id_vars=["group", "gene_id"])
        sailfish_melt["log10_value"] = np.log10(sailfish_melt["value"] + 1)
        
        feature_melt = pd.melt(feat_df,
                               id_vars=["group", "gene_id"])
        feature_melt["log10_value"] = np.log10(feature_melt["value"] + 1)

        feature_melt["variable"] = [fx.rstrip("_dedup") for fx in feature_melt["variable"]]
        melt_merge = pd.merge(sailfish_melt, feature_melt,
                              on=["variable", "gene_id"])
        melt_merge.columns = ["sailfish_group", "gene_id", "Sample", "sailfish_value",
                              "sailfish_log10_count", "feature_group", "feature_value",
                              "feature_log10_count"]

        return melt_merge


class ExpressionFacetGrid(object):
    '''
    Renderer for plotting FacetGrid of
    spike in transcript expression
    '''

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, data, path):

        ax = sns.FacetGrid(data, col="gene_id", col_wrap=5,
                           hue="gene_id")
        ax = (ax.map(plt.scatter, "sailfish_log10_count",
                     "feature_log10_count").set_titles("{col_name}"))
        ax.set(xlabel="Sailfish\nlog10 counts",
               ylabel="featureCounts\nlog10 counts")

        return ResultBlocks(ResultBlock(
            '''#$mpl %i$#\n''' % ax.figure,
            title="SpikeInPlot"))
