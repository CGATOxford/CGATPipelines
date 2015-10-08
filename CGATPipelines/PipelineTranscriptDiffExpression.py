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
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

from CGATPipelines.Pipeline import cluster_runnable

from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import rpy2.robjects as ro

import pandas as pd
import numpy as np


@cluster_runnable
def runSleuth(design, base_dir, model, contrast, outfile, counts, tpm, fdr):
    ''' run sleuth. Note: all samples in the design table must also
    have a directory with the same name in `base_dir` with kallisto
    results in a file called abundance.h5'''

    outfile_prefix = P.snip(outfile, ".tsv")

    Design = Expression.ExperimentalDesign(design)
    exp = Expression.DEExperiment_Sleuth()
    res = exp.run(base_dir, Design, model, contrast, outfile_prefix,
                  counts, tpm, fdr)

    res.getResults(fdr)

    res.plotMA(contrast, outfile_prefix)
    res.plotVolcano(contrast, outfile_prefix)

    res.table.to_csv(outfile, sep="\t", index=False)


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
        low_fraction = df[df['fraction_bin'] < 0.01].index.tolist()

        for transcript in low_fraction:
            outf.write("%s\t%s\n" % (transcript, "low_kmers"))

        # identify transcript with poor accuracy of quantification
        low_accuracy = df[[abs(x) > 0.5 for x in
                           df['log2diff']]].index.tolist()
        for transcript in low_accuracy:
            outf.write("%s\t%s\n" % (transcript, "poor_accuracy"))
