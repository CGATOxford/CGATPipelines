import numpy as np
import pandas as pd
import itertools
import collections
import random
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.preprocessing import scale as sklearn_scale
from sklearn.decomposition import PCA as sklearnPCA
from rpy2.robjects import r as R
import rpy2.robjects.pandas2ri as py2ri
from CGATReport.ResultBlock import ResultBlocks, ResultBlock
import seaborn
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import CGATPipelines.PipelineTracks as PipelineTracks
from CGATReport.Utils import PARAMS


DATABASE = PARAMS.get('readqc_backend', PARAMS.get('sql_backend', 'sqlite:///./csvdb'))


class RnaseqqcTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class SampleOverlap(RnaseqqcTracker):

    table = "sailfish_transcripts"

    def __call__(self, *args):

        # get list of unique IDs
        statement = ("""SELECT DISTINCT sample_id FROM '%(table)s'""")
        sample_list = self.getValues(statement)

        # create a DataFrame to output results
        df_range = range(1, len(sample_list)+1)
        result_df = pd.DataFrame(0, index=df_range, columns=df_range)

        # get all data at once
        statement = ("""SELECT sample_id, transcript_id, TPM
        FROM %(table)s
        WHERE transcript_id != 'Transcript'
        AND TPM >= 100""")

        working_df = self.getDataFrame(statement)

        # pairwise comparison of samples with common transcripts
        for samples in itertools.combinations_with_replacement(sample_list, 2):
            # get list of expressed transcripts for sample1
            transcripts_s1 = working_df[(working_df.sample_id == samples[0])]['transcript_id']
            # get list of expressed transcripts for sample2
            transcripts_s2 = working_df[(working_df.sample_id == samples[1])]['transcript_id']
            # compute intersection
            number_in_common = len(set(transcripts_s1) & set(transcripts_s2))
            # and update dataframe containing results
            result_df.iat[samples[0]-1, samples[1]-1] = number_in_common
            # this is a symmetrical matrix
            result_df.iat[samples[1]-1, samples[0]-1] = number_in_common

        return result_df


class SampleHeatmap(RnaseqqcTracker):

    table = "transcript_quantification"
    py2ri.activate()

    def getTracks(self, subset=None):
        return ("all")

    def getCorrelations(self, dataframe):
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
        corr_frame: pandas.Core.DataFrame
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

    def getFactors(self, dataframe, factor):
        '''
        Get factor/experimental design levels
        from table
        '''

        statement = ("SELECT factor_value, sample_id FROM factor "
                     "WHERE factor = '%(factor)s';" % locals())

        factor_df = pd.DataFrame.from_dict(self.getAll(statement))

        merged = pd.merge(dataframe, factor_df,
                          left_index=True, right_on="sample_id",
                          how='outer')
        return merged

    def __call__(self, track, slice=None):
        statement = ("SELECT sample_id,transcript_id,TPM from %(table)s "
                     "WHERE transcript_id != 'Transcript';")
        df = pd.DataFrame.from_dict(self.getAll(statement))

        mdf = self.getCorrelations(df)

        mdf.columns = set(df["sample_id"])
        mdf.index = set(df["sample_id"])

        all_df = self.getFactors(mdf, 'replicate')
        return all_df


class TranscriptQuantificationHeatmap(object):

    def __init__(self, *args, **kwargs):
        pass

    def getColorBar(self, data):

        # factors are the columns after the total
        # number of samples
        factors = data.iloc[:,data.shape[0]:]
        unique = set(factors.iloc[:,0])

        # select a random set of colours from the xkcd
        # palette
        random.seed(5648546)
        xkcd = random.sample(seaborn.xkcd_rgb.keys(),
                             len(unique)) 
        col_dict = dict(zip(unique, xkcd))
        cols = []
        for i in range(0, len(factors.index)):
            cols.append(seaborn.xkcd_rgb[col_dict[factors.iloc[i, 0]]])

        return cols

    def __call__(self, data, path):

        colorbar = self.getColorBar(data)
        n_samples = data.shape[0]
        data = data.iloc[:,:n_samples]
        print colorbar
        ax = seaborn.clustermap(data, 
                                row_colors=colorbar)

        return ResultBlocks(ResultBlock(
            '''#$mpl %i$#\n''' % ax.cax.figure.number,
            title='ClusterMapPlot'))


class sampleMDS(RnaseqqcTracker):
    # to add:
    # - ability to use rlog or variance stabalising transformatio
    # - ability to change filter threshold fo rlowly expressed transcripts
    # - JOIN with design table to get further aesthetics for plotting
    #   E.g treatment, replicate, etc

    table = "sailfish_transcripts"

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


class samplePCA(RnaseqqcTracker):
    '''
    Perform Principal component analysis on dataframe of
    expression values using sklearn PCA function. Takes expression
    dataframe, logs transforms data and scales variables to unit variance
    before performing PCA.
    '''

    # to add:
    # - ability to use rlog or variance stabalising transformation instead log2
    # - ability to change filter threshold for lowly expressed transcripts

    components = 10
    table = "sailfish_transcripts"

    def pca(self):

        # remove WHERE when table cleaned up to remove header rows
        statement = ("""SELECT transcript_id, TPM, sample_id FROM %s
        where transcript_id != 'Transcript' """ % self.table)

        # fetch data
        df = self.getDataFrame(statement)

        # put dataframe so row=genes, cols = samples, cells contain TPM
        pivot_df = df.pivot('transcript_id', 'sample_id')['TPM']

        # filter dataframe to get rid of genes where TPM == 0 across samples
        filtered_df = pivot_df[pivot_df.sum(axis=1) > 0]

        # add +1 to counts and log transform data.
        logdf = np.log(filtered_df + 0.1)

        # Scale dataframe so variance =1 across rows
        logscaled = sklearn_scale(logdf, axis=1)

        # turn array back to df and add transcript id back to index
        logscaled_df = pd.DataFrame(logscaled)
        logscaled_df.index = list(logdf.index)

        # Now do the PCA - can change n_components
        sklearn_pca = sklearnPCA(n_components=self.components)
        sklearn_pca.fit(logscaled_df)

        index = logdf.columns

        return sklearn_pca, index


class samplePCAprojections(samplePCA):
    '''
    Perform Principal component analysis on dataframe of
    expression values using sklearn PCA function. Takes expression
    dataframe, logs transforms data and scales variables to unit variance
    before performing PCA.

    Arguments
    ---------
    dataframe: pandas.Core.DataFrame
    a dataframe containing gene IDs, sample IDs
    and gene expression values

    Returns
    -------
    dataframe : pandas.Core.DataFrame
    a dataframe of first(PC1) and second (PC2) pricipal components
    in columns across samples, which are across the rows. '''

    # to add:
    # - ability to use rlog or variance stabalising transformation instead log2
    # - ability to change filter threshold for lowly expressed transcripts

    def __call__(self, track,  slice=None):

        sklearn_pca, index = self.pca()

        # these are the principle componets row 0 = PC1, 1 =PC2 etc
        PC_df = pd.DataFrame(sklearn_pca.components_)
        PC_df = PC_df.T
        PC_df.columns = ["PC%i" % x for x in range(1, self.components+1)]
        PC_df.index = index

        # This is what want for ploting bar graph
        # y = sklearn_pca.explained_variance_ratio_

        factor_statement = '''select * from factor'''

        # fetch factor data-THIS NEEDS TO BE ADJUSTED IF FACTORS table corrected
        factor_df = self.getDataFrame(factor_statement)
        factor_df.set_index("sample_id", drop=True, inplace=True)

        full_df = PC_df.join(factor_df)

        return collections.OrderedDict({x: full_df[full_df['factor'] == x] for
                                        x in set(full_df['factor'].tolist())})


class samplePCAvariance(samplePCA):
    '''
    Perform Principal component analysis on dataframe of
    expression values using sklearn PCA function. Takes expression
    dataframe, logs transforms data and scales variables to unit variance
    before performing PCA.

    Arguments
    ---------
    dataframe: pandas.Core.DataFrame
    a dataframe containing gene IDs, sample IDs
    and gene expression values

    Returns
    -------
    dataframe : pandas.Core.DataFrame
    a dataframe of first(PC1) and second (PC2) pricipal components
    in columns across samples, which are across the rows. '''
    # to add:
    # - ability to use rlog or variance stabalising transformation instead log2
    # - ability to change filter threshold for lowly expressed transcripts

    def __call__(self, track,  slice=None):

        sklearn_pca, index = self.pca()

        variance = sklearn_pca.explained_variance_ratio_

        final_df = pd.DataFrame({"variance": variance,
                                 "PC": range(1, self.components+1)})

        return final_df


class BiasFactors(RnaseqqcTracker):
    table = "bias_binned_means"

    def getTracks(self):
        d = self.get("SELECT DISTINCT bias_factor FROM %(table)s")
        return ["GC_Content", "length"]
        # return tuple([x[0] for x in d])

    def __call__(self, track, slice=None):
        statement = """
        SELECT bin, sample_id, value
        FROM %(table)s
        WHERE bias_factor = '%(track)s'
        AND variable = 'LogTPM'"""
        # fetch data
        df = self.getDataFrame(statement)
        df.set_index("sample_id", drop=False, inplace=True)

        factor_statement = '''select * from factors'''
        factor_df = self.getDataFrame(factor_statement)
        factor_df.set_index("sample_name", drop=True, inplace=True)
        factor_df.index.name = "sample_id"

        print factor_df.head()
        print df.head()

        full_df = df.join(factor_df)

        return full_df
        return collections.OrderedDict({x: full_df[full_df['factor'] == x] for
                                        x in set(full_df['factor'].tolist())})

        # TS: this should be replaced with a merge with the table of
        # experiment information
        # df2 = pd.DataFrame(map(lambda x: x.split("-"), df['sample']))
        # df2.columns = ["id_"+str(x) for x in range(1, len(df2.columns)+1)]

        # merged = pd.concat([df, df2], axis=1)
        # merged.index = ("all",)*len(merged.index)
        # merged.index.name = "track"


class ExpressionDistribution(RnaseqqcTracker):

    table = "sailfish_transcripts"

    def __call__(self, track, slice=None):
        statement = """SELECT sample_id, transcript_id, TPM
        FROM %(table)s WHERE transcript_id != 'Transcript'"""

        df = pd.DataFrame.from_dict(self.getAll(statement))
        c = 0.1
        df['logTPM'] = df['TPM'].apply(lambda x: np.log2(c + x))

        return df


class SampleOverlapsExpress(RnaseqqcTracker):
    '''
    Tracker class to compute overlap of expression for each
    sample on a pair-wise basis.  Returns a table of
    sample x sample overlaps, where the overlap is the
    number of common genes expressed in each pair of
    samples.
    '''

    table = "sailfish_transcripts"

    def __call__(self, track, slice=None):
        statement = """SELECT sample_id, transcript_id
        FROM %(table)s
        WHERE TPM >= 100;"""

        df = pd.DataFrame.from_dict(self.getAll(statement))

        overlaps = self.getOverlaps(dataframe=df)
        return overlaps

    def getOverlaps(self, dataframe):
        '''
        Pass in a dataframe of samples and
        expressed genes > threshold.
        Return an nxn dataframe of sample
        overlaps
        '''
        dataframe.index = dataframe["sample_id"]
        samples = set(dataframe.index)
        pairs = itertools.combinations_with_replacement(iterable=samples,
                                                        r=2)
        _df = pd.DataFrame(columns=samples, index=samples)
        _df.fillna(0.0, inplace=True)

        for comb in pairs:
            s1, s2 = comb
            s1_gene = set(dataframe.loc[s1]["transcript_id"])
            s2_gene = set(dataframe.loc[s2]["transcript_id"])
            gene_intersect = s1_gene.intersection(s2_gene)
            size = len(gene_intersect)
            _df.loc[s1, s2] = size
            _df.loc[s2, s1] = size

        return _df


class ThreePrimeBias(RnaseqqcTracker):
    '''
    Generates a dataframe of  mean read depth at each site in 3000bp from
    the 3' end.

    Arguments
    ---------
    threeprimebiasprofiles: str
    the name of an SQL database table containing mean read
    depth at these 3000 sites, plus upstream, downstream and intronic regions
    (these are not used)

    Returns
    -------
    df : pandas.Core.DataFrame
    a dataframe of showing bin (1 - 3000 with 3000 at the 3' end) and mean read
    count for only for the first 3000bp of each transcript'

    '''

    table = "threeprimebiasprofiles"

    def getTracks(self):
        d = self.get("""SELECT DISTINCT
                    track FROM threeprimebiasprofiles""")
        return tuple([x[0] for x in d])

    def __call__(self, track):
        statement = """
        SELECT bin, region, counts
        FROM %(table)s
        WHERE track = '%(track)s'
        AND region = 'exonsLast3000bp_zoomedTo3000bp'
        """
        df = self.getDataFrame(statement)
        df.bin -= 1000
        # reindexes bins as downstreem region not included
        return df

# class ExpressionDistributionNotR(RnaseqqcTracker, SingleTableTrackerColumns):
#    table = "sailfish_transcripts"
#    column = "transcript_id"
#    exclude_columns = "RPKM"

#    def __call__(self, track, slice=None):
#        statement = ("SELECT sample_id, transcript_id, RPKM FROM %(table)s WHERE transcript_id != 'Transcript'")
#        df = pd.DataFrame.from_dict(self.getAll(statement))
#        c = 0.0000001
#        df['log2rpkm'] = df['RPKM'].apply(lambda x: np.log2(c + x))
#        pivot = df.pivot(index='sample_id', columns='transcript_id', values='log2rpkm')

#        return pivot

# cgatreport-test -t ExpressionDistribution -r density-plot


class MappingTracker(TrackerSQL):
    """Base class for trackers from mapping report used for mapping context below"""

class MappingContext(MappingTracker, SingleTableTrackerRows):
    table = "context_stats"

class ProteinContext(RnaseqqcTracker):
    table = "context_stats"
    select ['protein_coding','non_stop_decay','nonsense_mediated_decay','IG_C_gene','IG_J_gene','IG_V_gene','TR_C_gene','TR_J_gene','TR_V_gene','processed_pseudogene','unprocessed_pseudogene','transcribed_processed_pseudogene','transcribed_unprocessed_pseudogene','translated_processed_pseudogene','unitary_pseudogene','IG_C_pseudogene','IG_J_pseudogene','IG_V_pseudogene','TR_J_pseudogene','TR_V_pseudogene','polymorphic_pseudogene','pseudogene']
    select=",".join(select)

    def getTracker(self, subset=None):
        retun ("all")

    def__call__(self, track, slice=None)
    statement = ("SELECT track,%(select)s FROM %(table)s")
    df = self.getDataFrame(statement)
    df.set_index("track", drop=True, inplace=True)
    df.index.name = "track"
    df = pd.melt(df, id_vars=['track'], var_name ='context')
    return df

class lncRNAContext(RnaseqqcTracker):
    table = "context_stats"
    select ['lincRNA','antisense','sense_overlapping','sense_intronic','_3prime_overlapping_ncrna','processed_transcript']
    select=",".join(select)

class ncRNAContext(RnaseqqcTracker):
    table = "context_stats"
    select ['snRNA','snoRNA','tRNA','Mt_rRNA', 'Mt_tRNA']
    select=",".join(select)

class riboRNAContext(RnaseqqcTracker):
    table = "context_stats"
    select ['rRNA','ribosomal_coding']
    select=",".join(select)



