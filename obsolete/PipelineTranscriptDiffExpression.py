'''
PipelineTranscriptDiffExpression.py - Utility functions for
pipeline_transcriptdiffexpression.py
==============================================================

:Author: Toms Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Code
----

'''

import CGAT.Expression as Expression
import CGAT.Counts as Counts
import CGAT.IOTools as IOTools

import CGATPipelines.Pipeline as P

from CGATPipelines.Pipeline import cluster_runnable

from rpy2.robjects import r as R

import pandas as pd
import numpy as np
import sqlite3
import os


def connect(database, annotations_database):
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(database)
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        annotations_database)
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


@cluster_runnable
def runSleuth(design, base_dir, model, contrasts, outfile, counts, tpm,
              fdr, lrt=False, reduced_model=None):
    ''' run sleuth. Note: all samples in the design table must also
    have a directory with the same name in `base_dir` with kallisto
    results in a file called abundance.h5'''

    outfile_prefix = P.snip(outfile, ".tsv")

    Design = Expression.ExperimentalDesign(design)
    exp = Expression.DEExperiment_Sleuth()

    res = exp.run(Design, base_dir, model, contrasts, outfile_prefix,
                  counts, tpm, fdr, lrt, reduced_model)

    res.getResults(fdr)
    for contrast in set(res.table['contrast']):
        res.plotMA(contrast, outfile_prefix)
        res.plotVolcano(contrast, outfile_prefix)

    res.table.to_csv(outfile, sep="\t", index=False)


@cluster_runnable
def runSleuthAll(samples, base_dir, counts, tpm):
    ''' run sleuth for all samples to obtain counts and tpm tables

    Note: all samples in the design table must also
    have a directory with the same name in `base_dir` with kallisto
    results in a file called abundance.h5
    '''

    design = pd.DataFrame({
        "group": ([0, 1] * ((len(samples) + 1) / 2))[0:len(samples)],
        "include": [1, ] * len(samples),
        "pair": [0, ] * len(samples)})

    design.index = samples

    Design = Expression.ExperimentalDesign(design)
    exp = Expression.DEExperiment_Sleuth()

    res = exp.run(Design, base_dir, counts=counts, tpm=tpm,
                  model="~group", dummy_run=True)


@cluster_runnable
def makeExpressionSummaryPlots(counts_inf, design_inf, logfile):
    ''' use the plotting methods for Counts object to make summary plots'''

    with IOTools.openFile(logfile, "w") as log:

        plot_prefix = P.snip(logfile, ".log")

        # need to manually read in data as index column is not the first column
        counts = Counts.Counts(pd.read_table(counts_inf, sep="\t"))
        counts.table.set_index(["transcript_id"])

        design = Expression.ExperimentalDesign(design_inf)

        # make certain counts table only include samples in design
        counts.restrict(design)

        cor_outfile = plot_prefix + "_pairwise_correlations.png"
        pca_var_outfile = plot_prefix + "_pca_variance.png"
        pca1_outfile = plot_prefix + "_pc1_pc2.png"
        pca2_outfile = plot_prefix + "_pc3_pc4.png"
        heatmap_outfile = plot_prefix + "_heatmap.png"

        counts_log10 = counts.log(base=10, pseudocount=0.1, inplace=False)

        counts_highExp = counts_log10.clone()
        counts_highExp.table['order'] = counts_highExp.table.apply(
            np.mean, axis=1)
        counts_highExp.table.sort(["order"], ascending=0, inplace=True)
        counts_highExp.table = counts_highExp.table.iloc[0:500, :]
        counts_highExp.table.drop("order", axis=1, inplace=True)

        log.write("plot correlations: %s\n" % cor_outfile)
        counts_log10.plotPairwiseCorrelations(cor_outfile, subset=1000)

        log.write("plot pc3,pc4: %s\n" % pca1_outfile)
        counts_log10.plotPCA(design,
                             pca_var_outfile, pca1_outfile,
                             x_axis="PC1", y_axis="PC2",
                             colour="group", shape="group")

        log.write("plot pc3,pc4: %s\n" % pca2_outfile)
        counts_log10.plotPCA(design,
                             pca_var_outfile, pca2_outfile,
                             x_axis="PC3", y_axis="PC4",
                             colour="group", shape="group")

        log.write("plot heatmap: %s\n" % heatmap_outfile)
        counts_highExp.heatmap(heatmap_outfile)


@cluster_runnable
def identifyLowConfidenceTranscripts(infile, outfile):
    ''' identify transcripts which cannot be confidently quantified in
    the simulation '''

    df = pd.read_table(infile, sep="\t", index_col=0)

    with IOTools.openFile(outfile, "w") as outf:

        outf.write("%s\t%s\n" % ("transcript_id", "reason"))

        # identify transcript with low fraction of kmers - these show
        # poorer correlation between ground truth and esimated counts
        low_fraction = df[df['fraction_bin'] < 0.03].index.tolist()

        for transcript in low_fraction:
            outf.write("%s\t%s\n" % (transcript, "low_kmers"))

        # identify transcript with poor accuracy of quantification
        low_accuracy = df[[abs(x) > 0.585 for x in
                           df['log2diff_tpm']]].index.tolist()

        for transcript in low_accuracy:
            outf.write("%s\t%s\n" % (transcript, "poor_accuracy"))


@cluster_runnable
def mergeAbundanceCounts(infile, outfile, counts):
    ''' merge the abundance and simulation counts files for
    each simulation '''

    df_abund = pd.read_table(infile, sep="\t", index_col=0)
    df_counts = pd.read_table(counts, sep="\t", index_col=0)
    df_abund.columns = [x if x != "tpm" else "est_tpm"
                        for x in df_abund.columns]

    df_merge = pd.merge(df_abund, df_counts, left_index=True, right_index=True)
    df_merge.index.name = "id"
    df_merge.to_csv(outfile, sep="\t")


@cluster_runnable
def calculateCorrelations(infiles, outfile, bin_step=1):
    ''' calculate correlation across simulation iterations per transcript'''

    abund, kmers = infiles

    df_abund = pd.read_table(abund, sep="\t", index_col=0)
    df_kmer = pd.read_table(kmers, sep="\t", index_col=0)

    # this is hacky, it's doing all against all correlations for the
    # two columns and subsetting
    df_agg_tpm = df_abund.groupby(level=0)[[
        "est_tpm", "tpm"]].corr().ix[0::2, 'tpm']

    # drop the "read_count" level, make into dataframe and rename column
    df_agg_tpm.index = df_agg_tpm.index.droplevel(1)
    df_agg_tpm = pd.DataFrame(df_agg_tpm)
    df_agg_tpm.columns = ["tpm_cor"]

    df_agg_count = df_abund.groupby(level=0)[[
        "est_counts", "read_count"]].corr().ix[0::2, 'read_count']

    # drop the "read_count" level, make into dataframe and rename column
    df_agg_count.index = df_agg_count.index.droplevel(1)
    df_agg_count = pd.DataFrame(df_agg_count)
    df_agg_count.columns = ["counts_cor"]

    # merge and bin the unique fraction values
    df_agg = pd.merge(df_agg_count, df_agg_tpm,
                      left_index=True, right_index=True)
    df_final = pd.merge(df_kmer, df_agg, left_index=True, right_index=True)
    df_final['fraction_bin'] = (
        np.digitize(df_final["fraction_unique"] * 100, bins=list(range(0, 100, bin_step)),
                    right=True)) / 100.0

    # Multiply bin number by step size to get the fraction for the bin
    df_final['fraction_bin'] = df_final['fraction_bin'] * bin_step

    df_abund_tpm_sum = df_abund.groupby(level=0)["est_tpm", "tpm"].sum()
    df_abund_count_sum = df_abund.groupby(level=0)[
        "est_counts", "read_count"].sum()
    df_abund_sum = pd.merge(df_abund_tpm_sum, df_abund_count_sum,
                            left_index=True, right_index=True)

    df_final = pd.merge(df_final, df_abund_sum,
                        left_index=True, right_index=True)

    df_final['log2diff_tpm'] = np.log2(df_final['est_tpm'] /
                                       df_final['tpm'])
    df_final['log2diff_tpm_thres'] = [x if abs(x) < 2 else 2 * x / abs(x)
                                      for x in df_final['log2diff_tpm']]

    df_final['log2diff_counts'] = np.log2(df_final['est_counts'] /
                                          df_final['read_count'])
    df_final['log2diff_counts_thres'] = [x if abs(x) < 1 else x / abs(x)
                                         for x in df_final['log2diff_counts']]

    df_final.to_csv(outfile, sep="\t", index=True)


def loadSleuthTable(infile, outfile, transcript_info, gene_biotypes,
                    database, annotations_database):

    tmpfile = P.getTempFilename("/ifs/scratch/")

    table = os.path.basename(transcript_info)

    if gene_biotypes:
        where_cmd = "WHERE " + " OR ".join(
            ["gene_biotype = '%s'" % x
             for x in gene_biotypes.split(",")])
    else:
        where_cmd = ""

    select = """SELECT DISTINCT
        transcript_id, transcript_biotype, gene_id, gene_name
        FROM annotations.%(table)s
        %(where_cmd)s""" % locals()

    df1 = pd.read_table(infile, sep="\t")
    df1.set_index("transcript_id", drop=True, inplace=True)

    df2 = pd.read_sql(select, connect(database, annotations_database))
    df2.set_index("transcript_id", drop=False, inplace=True)

    df = df1.join(df2)
    df.to_csv(tmpfile, sep="\t", index=True)

    options = "--add-index=transcript_id"
    P.load(tmpfile, outfile, options=options)
    os.unlink(tmpfile)


def loadSleuthTableGenes(infile, outfile, gene_info, gene_biotypes,
                         database, annotations_database):

    tmpfile = P.getTempFilename("/ifs/scratch/")

    table = os.path.basename(gene_info)

    if gene_biotypes:
        where_cmd = "WHERE " + " OR ".join(
            ["gene_biotype = '%s'" % x
             for x in gene_biotypes.split(",")])
    else:
        where_cmd = ""

    select = """SELECT DISTINCT
        gene_id, gene_name
        FROM annotations.%(table)s
        %(where_cmd)s""" % locals()

    df1 = pd.read_table(infile, sep="\t")
    df1.set_index("test_id", drop=False, inplace=True)

    df2 = pd.read_sql(select, connect(database, annotations_database))
    df2.set_index("gene_id", drop=False, inplace=True)

    df = df1.join(df2)
    df.to_csv(tmpfile, sep="\t", index=True)

    options = "--add-index=gene_id"
    P.load(tmpfile, outfile, options=options)
    os.unlink(tmpfile)


def convertFromFish(infile, outfile):
    ''' convert sailfish/salmon output to Sleuth compatible h5 file'''

    infile = os.path.dirname(infile)

    convert = R('''library(wasabi)
    prepare_fish_for_sleuth("%(infile)s", force=TRUE)
    ''' % locals())

    convert
