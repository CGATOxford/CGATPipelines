import numpy as np
import pandas as pd
import itertools
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
from CGATReport.Utils import PARAMS
import CGATPipelines.Pipeline as Pipeline
import matplotlib.pyplot as plt

DATABASE = PARAMS.get('readqc_backend',
                      PARAMS.get('sql_backend', 'sqlite:///./csvdb'))


class imagesTracker(TrackerImages):

    '''Convience Tracker for globbing images for gallery plot'''
    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]


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
        df_range = list(range(1, len(sample_list)+1))
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

        return result_df.reset_index().set_index("sample_id")


class SampleHeatmap(RnaseqqcTracker):

    table = "sailfish_transcripts"
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
        pivot = dataframe.pivot(index="sample_name",
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

    def getFactors(self, dataframe):
        '''Get factor/experimental design levels from table
        '''

        statement = ("SELECT factor_value, sample_name, factor "
                     "FROM factors AS f "
                     "JOIN samples AS s "
                     "ON f.sample_id = s.id "
                     "WHERE factor != 'genome'" % locals())

        factor_df = self.getDataFrame(statement)
        merged = pd.merge(dataframe,
                          factor_df,
                          left_index=True,
                          right_on="sample_name",
                          how='outer')
        return merged

    def __call__(self, track, slice=None):

        statement = ("SELECT s.sample_name, t.transcript_id, t.TPM "
                     "FROM %(table)s AS t, samples AS s "
                     "WHERE transcript_id != 'Transcript' "
                     "AND t.sample_id = s.id")
        df = self.getDataFrame(statement)
        mdf = self.getCorrelations(df)

        mdf.columns = set(df["sample_name"])
        mdf.index = set(df["sample_name"])

        all_df = self.getFactors(mdf)
        return all_df.set_index("factor")


class TranscriptQuantificationHeatmap(object):

    def __init__(self, *args, **kwargs):
        pass

    def getColorBar(self, data):

        # factors are the columns after the total number of samples
        factors = data.iloc[:, data.shape[0]:]
        unique = set(factors.iloc[:, 0])

        # select a random set of colours from the xkcd palette
        random.seed(5648546)
        xkcd = random.sample(list(seaborn.xkcd_rgb.keys()),
                             len(unique))
        col_dict = dict(list(zip(unique, xkcd)))
        cols = []
        for i in range(0, len(factors.index)):
            cols.append(seaborn.xkcd_rgb[col_dict[factors.iloc[i, 0]]])

        return cols, factors, unique, xkcd

    def __call__(self, data, path):

        colorbar, factors, unique, xkcd = self.getColorBar(data)
        n_samples = data.shape[0]
        data = data.iloc[:, :n_samples]
        col_dict = dict(list(zip(unique, xkcd)))

        print(data.head())
        seaborn.set(font_scale=.5)
        ax = seaborn.clustermap(data,
                                row_colors=colorbar, col_colors=colorbar)
        plt.setp(ax.ax_heatmap.yaxis.set_visible(False))

        for label in unique:
            ax.ax_col_dendrogram.bar(
                0, 0, color=seaborn.xkcd_rgb[col_dict[label]],
                label=label, linewidth=0)
        ax.ax_col_dendrogram.legend(loc="center", ncol=len(unique))

        return ResultBlocks(ResultBlock(
            '''#$mpl %i$#\n''' % ax.cax.figure.number,
            title='ClusterMapPlot'))


class sampleMDS(RnaseqqcTracker):

    def __call__(self, track,  slice=None):

        # remove WHERE when table cleaned up to remove header rows
        statement = (
            "SELECT transcript_id, TPM, sample_id FROM sailfish_transcripts")

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

        factors_df = self.getDataFrame(
            "SELECT * FROM factors WHERE factor != 'genome'")

        merged_df = pd.merge(pos, factors_df,
                             left_on="sample", right_on="sample_id")
        return merged_df.reset_index().set_index("factor")


class SamplePCA(RnaseqqcTracker):
    '''
    Perform Principal component analysis on dataframe of
    expression values using sklearn PCA function. Takes expression
    dataframe, logs transforms data and scales variables to unit variance
    before performing PCA.
    '''

    # to add:
    # - ability to use rlog or variance stabalising transformation instead log2
    # - ability to change filter threshold for lowly expressed transcripts

    n_components = None
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
        sklearn_pca = sklearnPCA(n_components=self.n_components)
        sklearn_pca.fit(logscaled_df)

        index = logdf.columns

        return sklearn_pca, index


class SamplePCAProjections(SamplePCA):
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

        PC_df.columns = ["PC%i" % x for x in
                         range(1, sklearn_pca.n_components_ + 1)]
        PC_df.index = index

        # This is what want for ploting bar graph
        # y = sklearn_pca.explained_variance_ratio_

        factor_df = self.getDataFrame(
            "SELECT * FROM factors where factor != 'genome'")
        factor_df.set_index("sample_id", inplace=True)
        if "replicate" in set(factor_df.factor):
            replicate_df = factor_df[factor_df.factor == "replicate"]
            replicate_df.columns = ["x", "replicate"]
            factor_df = factor_df.join(replicate_df).drop("x", axis=1)
            factor_df = factor_df[factor_df.factor != "replicate"]
        else:
            factor_df["replicate"] = "R0"

        # remove factors with a single level
        full_df = PC_df.join(factor_df)
        return full_df.reset_index().set_index(["factor", "sample_id"])


class SamplePCAVariance(SamplePCA):
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

        final_df = pd.DataFrame(
            {"variance": variance,
             "PC": list(range(1, sklearn_pca.n_components_ + 1))})

        return final_df


class BiasFactors(RnaseqqcTracker):

    bias_factor = ""

    def __call__(self, track, slice=None):

        statement = """
        SELECT bin, sample_id, value
        FROM bias_binned_means
        WHERE variable = 'LogTPM_norm' AND bias_factor='%(bias_factor)s'
        """

        # fetch data
        df = self.getDataFrame(statement)
        df.set_index("sample_id", drop=False, inplace=True)

        factor_statement = '''select * from factors where factor != "genome"'''
        factor_df = self.getDataFrame(factor_statement)

        full_df = pd.merge(df, factor_df, left_on="sample_id",
                           right_on="sample_id")

        full_df = full_df.reset_index().set_index("factor")
        full_df.index.names = ["track"]

        return full_df


class BiasFactorsGCContent(BiasFactors):
    bias_factor = "GC_Content"


class BiasFactorsLength(BiasFactors):
    bias_factor = "length"


class BiasFactorsTT(BiasFactors):
    bias_factor = "TT"


class BiasFactorsAA(BiasFactors):
    bias_factor = "AA"


class BiasFactorsCC(BiasFactors):
    bias_factor = "CC"


class BiasFactorsGG(BiasFactors):
    bias_factor = "GG"


class ExpressionDistribution(RnaseqqcTracker):

    def __call__(self, track, slice=None):
        statement = """
        SELECT sample_id, transcript_id, TPM
        FROM sailfish_transcripts"""

        df = self.getDataFrame(statement)

        df['logTPM'] = df['TPM'].apply(lambda x: np.log2(x + 0.1))

        factor_statement = '''
        select *
        FROM factors
        JOIN samples
        ON factors.sample_id = samples.id
        WHERE factor != "genome"'''

        factor_df = self.getDataFrame(factor_statement)

        full_df = pd.merge(df, factor_df, left_on="sample_id",
                           right_on="sample_id")

        print(full_df.reset_index().set_index("factor").shape)
        # trying to find the point at which the error occurs!
        return full_df.reset_index().set_index("factor").tail(10355125)


class SampleOverlapsExpress(RnaseqqcTracker):
    '''
    Tracker class to compute overlap of expression for each
    sample on a pair-wise basis.  Returns a table of
    sample x sample overlaps, where the overlap is the
    number of common genes expressed in each pair of
    samples.
    '''

    def __call__(self, track, slice=None):
        statement = """SELECT sample_name, transcript_id
        FROM sailfish_transcripts
        JOIN samples on sample_id = id
        WHERE TPM >= 100;"""

        df = pd.DataFrame.from_dict(self.getAll(statement))

        overlaps = self.getOverlaps(df)
        return overlaps

    def getOverlaps(self, dataframe):
        '''
        Pass in a dataframe of samples and
        expressed genes > threshold.
        Return an nxn dataframe of sample
        overlaps
        '''
        dataframe.index = dataframe["sample_name"]
        samples = set(dataframe.index)
        pairs = itertools.combinations_with_replacement(iterable=samples,
                                                        r=2)
        df = pd.DataFrame(columns=samples, index=samples)
        df.fillna(0.0, inplace=True)

        for comb in pairs:
            s1, s2 = comb
            s1_gene = set(dataframe.loc[s1]["transcript_id"])
            s2_gene = set(dataframe.loc[s2]["transcript_id"])
            gene_intersect = s1_gene.intersection(s2_gene)
            size = len(gene_intersect)
            df.loc[s1, s2] = size
            df.loc[s2, s1] = size

        print(df.head())
        return df


class ThreePrimeBias(RnaseqqcTracker):
    '''
    Generates a dataframe of mean read depth at each site in 3000bp from
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


class MappingTracker(TrackerSQL):
    """Base class for trackers from mapping report"""


class MappingContext(MappingTracker, SingleTableTrackerRows):
    table = "context_stats"


class MappingContentRNA(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        statement = "SELECT rRNA, total, track FROM context_stats"
        df = self.getDataFrame(statement)
        df.set_index("track", drop=True, inplace=True)
        df.index.name = "track"

        def norm_row(array):
            total = array['total']
            return [float(x)/total for x in array]

        df = df.apply(axis=1, func=norm_row)
        df['other'] = df['total'] - df['rRNA']
        df.drop('total', axis=1, inplace=True)
        return df


class MappingAltContent(RnaseqqcTracker):
    table = "altcontext_stats"

    def __call__(self, track, slice=None):

        statement = "SELECT * FROM %(table)s"
        df = self.getDataFrame(statement)
        df.set_index("track", drop=True, inplace=True)
        df.index.name = "track"

        def norm_row(array):
            total = array['total']
            return [float(x)/total for x in array]

        df = df.apply(axis=1, func=norm_row)
        other_columns = [x for x in df.columns if "total" not in x]
        df['other'] = df['total'] - df.loc[:, other_columns].sum(axis=1)
        df.drop('total', axis=1, inplace=True)
        return df


class ProteinContext(RnaseqqcTracker):
    table = "context_stats"
    select = ['protein_coding', 'non_stop_decay', 'nonsense_mediated_decay',
              'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'TR_C_gene', 'TR_J_gene',
              'TR_V_gene', 'processed_pseudogene', 'unprocessed_pseudogene',
              'transcribed_processed_pseudogene',
              'transcribed_unprocessed_pseudogene',
              'translated_processed_pseudogene', 'unitary_pseudogene',
              'IG_C_pseudogene',
              'IG_J_pseudogene', 'IG_V_pseudogene', 'TR_J_pseudogene',
              'TR_V_pseudogene',
              'polymorphic_pseudogene', 'pseudogene']

    def getTracker(self, subset=None):
        return ("all")

    def __call__(self, track, slice=None):

        select = ",".join(self.select)

        statement = "SELECT track, %(select)s FROM %(table)s"
        df = self.getDataFrame(statement)
        df.set_index("track", drop=True, inplace=True)
        df.index.name = "track"
        df = pd.melt(df, id_vars=['track'], var_name='context')
        return df


class lncRNAContext(ProteinContext):
    select = ['lincRNA', 'antisense', 'sense_overlapping', 'sense_intronic',
              '_3prime_overlapping_ncrna', 'processed_transcript']


class ncRNAContext(ProteinContext):
    select = ['snRNA', 'snoRNA', 'tRNA', 'Mt_rRNA', 'Mt_tRNA']


class riboRNAContext(ProteinContext):
    select = ['rRNA', 'ribosomal_coding']


class TopGenes(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        exp_statement = """
        SELECT TPM, gene_id, sample_name
        FROM sailfish_genes AS A
        JOIN samples AS B
        ON A.sample_id = B.id"""

        exp_df = self.getDataFrame(exp_statement)

        factors_statement = '''
        SELECT factor, factor_value, sample_name
        FROM samples AS A
        JOIN factors AS B
        ON A.id = B.sample_id
        WHERE factor != 'genome'
        '''

        factors_df = self.getDataFrame(factors_statement)

        exp_df['TPM'] = exp_df['TPM'].astype(float)
        exp_df_pivot = pd.pivot_table(exp_df, values=["TPM"],
                                      index="gene_id", columns="sample_name")

        rowSums = exp_df_pivot.apply(axis=1, func=sum)
        top_genes_df = exp_df_pivot.ix[
            rowSums.sort_values(ascending=False)[0:20].index]

        def z_norm(array):
            mu = np.mean(array)
            sd = np.std(array)
            return [(x-mu)/sd for x in array]

        z_norm_df = top_genes_df.apply(axis=0, func=z_norm)
        z_norm_df['gene_id'] = z_norm_df.index

        z_norm_melted_df = pd.melt(z_norm_df, id_vars="gene_id")
        z_norm_melted_df.set_index('sample_name', inplace=True, drop=False)

        merged_df = pd.merge(z_norm_melted_df, factors_df,
                             left_on="sample_name", right_on="sample_name")

        return merged_df.reset_index().set_index("factor")


class GenesOfInterest(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        exp_statement = """
        SELECT TPM, gene_id, sample_name
        FROM sailfish_genes AS A
        JOIN samples AS B
        ON A.sample_id = B.id"""

        exp_df = self.getDataFrame(exp_statement)

        factors_statement = '''
        SELECT factor, factor_value, sample_name
        FROM samples AS A
        JOIN factors AS B
        ON A.id = B.sample_id
        WHERE factor != 'genome'
        '''

        factors_df = self.getDataFrame(factors_statement)

        merged_df = pd.merge(exp_df, factors_df,
                             left_on="sample_name", right_on="sample_name")

        genes = Pipeline.asList(Pipeline.peekParameters(
            ".", "pipeline_rnaseqqc.py")['genes_of_interest'])

        interest_df = merged_df[merged_df['gene_id'].isin(genes)]

        interest_df['TPM'] = interest_df['TPM'].astype(float)

        return interest_df.reset_index().set_index("factor")


class PicardThreePrimeBias(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        statement = '''
        SELECT MEDIAN_3PRIME_BIAS, MEDIAN_5PRIME_BIAS,
        MEDIAN_5PRIME_TO_3PRIME_BIAS, track
        FROM picard_rna_metrics'''

        df = self.getDataFrame(statement)

        # remove the ".hisat" suffix from track
        df['track'] = [x.split(".")[0] for x in df['track']]

        factors_statement = '''
        SELECT factor, factor_value, sample_name
        FROM samples AS A
        JOIN factors AS B
        ON A.id = B.sample_id
        WHERE factor != 'genome'
        '''

        factors_df = self.getDataFrame(factors_statement)

        merged_df = pd.merge(df, factors_df,
                             left_on="track", right_on="sample_name")

        final_df = pd.melt(merged_df, id_vars=[
            "track", "factor", "factor_value", "sample_name"])

        return final_df.set_index('factor')


class PicardStrandBias(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        statement = '''
        SELECT CORRECT_STRAND_READS, INCORRECT_STRAND_READS, track
        FROM picard_rna_metrics'''

        df = self.getDataFrame(statement)

        # remove the ".hisat" suffix from track
        df['track'] = [x.split(".")[0] for x in df['track']]
        df['correct_strand_perc'] = 100 * (
            df['CORRECT_STRAND_READS'] / (
                df['CORRECT_STRAND_READS'] + df['INCORRECT_STRAND_READS']))

        return df


class Factors(RnaseqqcTracker):

    def __call__(self, track, slice=None):

        select = '''SELECT factor, factor_value, sample_name
        FROM factors JOIN samples ON sample_id=id'''

        df = self.getDataFrame(select)
        pivot_df = df.pivot(
            index="sample_name", values="factor_value", columns="factor")

        return pivot_df
