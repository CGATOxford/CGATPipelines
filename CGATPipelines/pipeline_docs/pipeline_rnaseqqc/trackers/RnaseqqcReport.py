import re
import glob
import pandas as pd
import numpy as np
from sklearn import manifold
from sklearn.metrics import euclidean_distances
import rpy2.robjects as ro
from rpy2.robjects import r as R
import rpy2.robjects.pandas2ri as py2ri
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import CGATPipelines.PipelineTracks as PipelineTracks

###################################################################
###################################################################
# parameterization

EXPORTDIR = P.get('readqc_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('readqc_datadir', P.get('datadir', '.'))
DATABASE = P.get('readqc_backend', P.get('sql_backend', 'sqlite:///./csvdb'))

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################


TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.sra" % DATADIR), "(\S+).sra") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % DATADIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % DATADIR), "(\S+).fastq.1.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("*.csfasta.gz"), "(\S+).csfasta.gz")

###########################################################################


class RnaseqqcTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

##############################################################
##############################################################
##############################################################


class SampleHeatmap(RnaseqqcTracker):
    table = "transcript_quantification"
    py2ri.activate()

    def getTracks(self, subset=None):
        return ("all")

    def getCurrentRDevice(self):

        '''return the numerical device id of the
        current device'''

        return R["dev.cur"]()[0]

    def hierarchicalClustering(self, dataframe):
        '''
        Perform hierarchical clustering on a
        dataframe of expression values

        Arguments
        ---------
        dataframe: pandas.Core.DataFrame
          a dataframe containing gene IDs, sample IDs
          and gene expression values

        Returns
        -------
        correlations: pandas.Core.DataFrame
          a dataframe of a pair-wise correlation matrix
          across samples.  Uses the Pearson correlation.
        '''

        # set sample_id to index
        pivot = dataframe.pivot(index="sample_id",
                                columns="transcript_id",
                                values="TPM")
        transpose = pivot.T
        # why do I have to resort to R????
        r_df = py2ri.py2ri_pandasdataframe(transpose)
        R.assign("p.df", r_df)
        R('''p.mat <- apply(p.df, 2, as.numeric)''')
        R('''cor.df <- cor(p.mat)''')
        r_cor = R["cor.df"]
        py_cor = py2ri.ri2py_dataframe(r_cor)
        corr_frame = py_cor

        return corr_frame

    def __call__(self, track, slice=None):
        statement = ("SELECT sample_id,transcript_id,TPM from %(table)s "
                     "WHERE transcript_id != 'Transcript';")
        df = pd.DataFrame.from_dict(self.getAll(statement))
        # insert clustering function here

        mdf = self.hierarchicalClustering(df)
        mdf.columns = set(df["sample_id"])
        mdf.index = set(df["sample_id"])
        r_cor = py2ri.py2ri_pandasdataframe(mdf)
        R.assign("cor.mat", r_cor)

        R.x11()
        R('''suppressPackageStartupMessages(library(gplots))''')
        R('''suppressPackageStartupMessages(library(RColorBrewer))''')
        R('''hmcol <- colorRampPalette(c("#FFFF00", "#7A378B"))''')
        R('''heatmap.2(as.matrix(cor.mat), trace="none",'''
          '''col=hmcol)''')

        return odict((("Sum absolute covariance",
                       "#$rpl %i$#" % self.getCurrentRDevice()),))


class sampleMDS(RnaseqqcTracker):
    # to add:
    # - parameterise dissimilarity so we can plot
    #    euclidean & 1-cor(spearman's?)
    # - JOIN with design table to get further aesthetics for plotting
    #   E.g treatment, replicate, etc

    table = "transcript_quantification"

    def __call__(self, track,  slice=None):

        # remove WHERE when table cleaned up to remove header rows
        statement = (
            "SELECT transcript_id, TPM, sample_id FROM %(table)s "
            "where transcript_id != 'Transcript'")

        # fetch data
        df = pd.DataFrame.from_dict(self.getAll(statement))

        df = df.pivot('transcript_id', 'sample_id')['TPM']

        # calculate dissimilarities
        similarities = euclidean_distances(df.transpose())

        # run MDS
        mds = manifold.MDS(n_components=2, max_iter=3000,
                           eps=1e-9, dissimilarity="precomputed", n_jobs=1)
        mds = mds.fit(similarities)
        pos = pd.DataFrame(mds.embedding_)

        pos.columns = ["MD1", "MD2"]
        pos['sample'] = df.columns

        return pos


class CorrelationSummaryA(RnaseqqcTracker):
    table = "binned_means_correlation"
    select = ["AA", "AT", "AC", "AG"]
    select = ",".join(select)

    def __call__(self, track,  slice=None):
        statement = ("SELECT sample,%(select)s FROM %(table)s")
        # fetch data
        df = pd.DataFrame.from_dict(self.getAll(statement))
        df['sample'] = [x.replace("_quant.sf", "") for x in df['sample']]
        df = pd.melt(df, id_vars="sample")
        df2 = pd.DataFrame(map(lambda x: x.split("-"), df['sample']))
        df2.columns = ["id_"+str(x) for x in range(1, len(df2.columns)+1)]
        merged = pd.concat([df, df2], axis=1)
        merged.index = ("all",)*len(merged.index)
        merged.index.name = "track"
        return merged


class GradientSummaryA(CorrelationSummaryA):
    table = "binned_means_gradients"


class CorrelationSummaryT(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["TA", "TT", "TC", "TG"]
    select = ",".join(select)


class GradientSummaryT(CorrelationSummaryT):
    table = "binned_means_gradients"


class CorrelationSummaryC(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["CA", "CT", "CC", "CG"]
    select = ",".join(select)


class GradientSummaryC(CorrelationSummaryC):
    table = "binned_means_gradients"


class CorrelationSummaryG(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["GA", "GT", "GC", "GG"]
    select = ",".join(select)


class GradientSummaryG(CorrelationSummaryG):
    table = "binned_means_gradients"


class CorrelationSummaryGC(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["GC_Content", "length"]
    select = ",".join(select)


class GradientSummaryGC(CorrelationSummaryGC):
    table = "binned_means_gradients"


class BiasFactors(RnaseqqcTracker):
    table = "binned_means"

    def getTracks(self):
        d = self.get("SELECT DISTINCT factor FROM %(table)s")
        return tuple([x[0] for x in d])

    def __call__(self, track, slice=None):
        statement = "SELECT * FROM %(table)s WHERE factor = '%(track)s'"
        # fetch data
        df = self.getDataFrame(statement)

        # TS: this should be replaces with a merge with the table of
        # experiment information
        df2 = pd.DataFrame(map(lambda x: x.split("-"), df['sample']))
        df2.columns = ["id_"+str(x) for x in range(1, len(df2.columns)+1)]

        merged = pd.concat([df, df2], axis=1)
        merged.index = ("all",)*len(merged.index)
        merged.index.name = "track"

        return merged
