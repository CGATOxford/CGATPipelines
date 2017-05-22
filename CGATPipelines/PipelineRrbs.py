'''
PipelineRrbs.py - Utility functions for analysing rrbs sequence data
=====================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra` files)
of different data (single end, paired end) with different mapping
algorithms (bowtie, tophat, stampy). This module provides utility
functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra` archive
might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

The module currently is able to deal with:

   * tophat mapping against genome
   * bowtie mapping against transcriptome, genome and junctions
   * bwa against genome
   * stampy against genome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''
import os
import re
import sqlite3
import collections
import itertools
import copy
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator
import CGAT.IndexedFasta as IndexedFasta
import pandas as pd
from CGATPipelines.Pipeline import cluster_runnable
import numpy as np
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects import pandas2ri
# AH: this causes an import error:
# AssertionError: PipelineRrbs scripts/modules - ImportError: Conflict when converting R symbol in the package "stats" to a Python symbol (format.perc -> format_perc while there is already format_perc)
# best only to import when needed in individual methods
# stats = importr('stats')


def hjoin(items):
    return "-".join(items)


def nullConvert(start, end, strand, seq):
    ''' define a convert function so the positions for sequences
    on the '-' strand remain relative to the start of the
    contig. Otherwise PysamIndexedFasta will perform this
    conversion(!)'''
    return start, end


def getCpGs(seq):
    '''return postions of CpGs'''
    seq = seq.upper()
    CpGs = [x.start() for x in re.finditer("CG", seq)]
    return CpGs


@cluster_runnable
def plotReadBias(infile, outfile):
    '''take the m-bias outfile from bismark and plot read position biases'''

    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    def dfFromTsvTable(lines, line, context):
        row = 2
        rows = []
        row_values = lines[line + row].strip().split("\t")
        header = [re.sub(" ", "_", re.sub("%", "perc", x))
                  for x in row_values]
        row += 1
        while lines[line + row].strip():
            print(lines[line + row].strip())
            row_values = lines[line + row].strip().split("\t")
            row_values = list(map(num, row_values))
            rows.append(row_values)
            row += 1
        df = pd.DataFrame(rows, columns=header)
        df['context'] = context
        return(df)

    with open(infile, "r") as f:
        lines = f.readlines()
        line = 0
        while line < len(lines):
            if lines[line].strip() == "CpG context":
                df = dfFromTsvTable(lines, line, "CpG")
                print(df)
            elif lines[line].strip() == "CHG context":
                df2 = dfFromTsvTable(lines, line, "CHG")
                print(df2)
            elif lines[line].strip() == "CHH context":
                df3 = dfFromTsvTable(lines, line, "CHH")
                print(df3)
            line += 1
    final_df = pd.concat([df, df2, df3])
    final_df.index = list(range(len(final_df['context'])))
    r_dataframe = pandas2ri.py2ri(final_df)

    plot1_out = (P.snip(outfile, ".read_position.tsv") +
                 "_read_position_methylation_bias.png")
    plot2_out = re.sub("_methylation_bias", "_count", plot1_out)

    title = re.sub(".fastq.*", "", os.path.basename(infile))
    plotMbiasReadPositionSummary = r("""
    library(ggplot2)
    library(reshape)
    require(scales)
    function(df){
    read_length = max(df$position)
    p = ggplot(df, aes(x=position, y=perc_methylation, col=context)) +
    geom_line() + geom_point(size=1) +
    facet_grid(context ~ ., scales="free") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5)) +
    ggtitle("%(title)s") +
    scale_colour_discrete(guide = FALSE) +
    xlab("Read position") + ylab("Methylation (percentage)")
    ggsave('%(plot1_out)s', plot = p, width = 5, height = 5)

    p2 = ggplot(df, aes(x=position,
    y=count_methylated+count_unmethylated, col=context)) +
    geom_line() + geom_point(size=1) +
    facet_grid(context ~ ., scales="free") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5)) +
    ggtitle("%(title)s") +
    scale_colour_discrete(guide = FALSE) +
    scale_y_continuous(labels = comma) +
    xlab("Read position") + ylab("Number of cytosines")
    ggsave('%(plot2_out)s', plot = p2, width = 5, height = 5)

    write.table(df, file="%(outfile)s", row.names = F, sep="\t")
    }""" % locals())

    plotMbiasReadPositionSummary(r_dataframe)


@cluster_runnable
def pandasMerge(infile1, infile2, outfile, merge_type, left, right,
                delim="\t", com="#"):
    '''function to merge two tsv files using pandas
    left and right are the columns to merge on'''

    def columnsPeak(infile):
        '''peak at columns to identify data type'''
        top = pd.read_csv(infile, sep=delim, comment=com, nrows=1)
        columns = top.columns.tolist()

        # strand, contig and cpgi are all strings, hence object
        # all other columns encoded as floats
        dtype_dict = {}

        for column in columns:
            if ((column.endswith("-unmeth") or
                 column.endswith("-meth") or
                 column.endswith("-perc") or
                 column == "read_position")):
                dtype_dict[column] = float
            elif (column == "strand" or
                  column == "contig" or
                  column == "cpgi"):
                dtype_dict[column] = object
            elif column == "position":
                dtype_dict[column] = float

        return dtype_dict

    def pandasRead(infile, peak=False):
        if peak:
            col_dtype = columnsPeak(infile)
            print(col_dtype)
            return pd.read_csv(infile, sep=delim, comment=com,
                               dtype=col_dtype, na_values="NA")
        else:
            return pd.read_csv(infile, sep=delim, comment=com)

    df1 = pandasRead(infile1, peak=True)
    df2 = pandasRead(infile2, peak=True)
    df1.set_index(left, drop=False, append=False, inplace=True)
    df2.set_index(right, drop=True, append=False, inplace=True)

    merged = df1.merge(df2, how=merge_type, left_index=True,
                       right_index=True, sort=False)
    merged.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def fasta2CpG(infile, outfile):
    '''perform an in-silico digest at MspI sites and return all CpGs
    thorughout the genome, whether they are in a MspI fragment
    and if so, what their read position is'''
    # to do: paramterise (digestion site, read length, PE/SE)

    # AH: Use FastaIterator
    fasta = FastaIterator.FastaIterator(IOTools.openFile(infile, "r"))
    temp_contig = next(fasta)

    outfile = IOTools.openFile(outfile, "w")
    outfile.write("contig\tposition\tstrand\tread_position\n")

    while len(temp_contig.seq) > 1:
        contig, contig_seq = temp_contig.title, temp_contig.sequence
        # find MspI sites
        iter_pos = re.finditer("[cC][cC][gG][gG]", contig_seq)
        pos_start = [x.start(0) for x in iter_pos]
        # loop through pairs of MspI sites and add positions to dictionary
        MspI_cpg_dict = {}
        for n in range(0, len(pos_start) - 1):
            frag_length = pos_start[n + 1] - pos_start[n]
            if frag_length > 25 and frag_length < 300:
                if frag_length >= 51:
                    # need to include additional base at the end of the read
                    # to check for final CG
                    read_f = contig_seq[pos_start[n] + 1:pos_start[n] + 53]
                    read_r = contig_seq[
                        pos_start[n + 1] - 49:pos_start[n + 1] + 3]
                else:
                    read_f = contig_seq[pos_start[n] + 1:pos_start[n + 1] + 1]
                    read_r = contig_seq[pos_start[n] + 3:pos_start[n + 1] + 3]
                # find positions of CGs in in silico reads
                f_cgs = re.finditer("[cC][gG]", read_f)
                r_cgs = re.finditer("[cC][gG]", read_r)
                for x in f_cgs:
                    # 1-based read position
                    read_pos = x.start(0) + 1
                    # contig position 1-based
                    contig_position = pos_start[n] + 1 + read_pos
                    MspI_cpg_dict[contig_position] = read_pos
                for x in r_cgs:
                    # 1-based read position
                    read_pos = len(read_r) - (x.start(0) + 1)
                    # contig position 1-based
                    contig_position = pos_start[n + 1] + 4 - read_pos
                    MspI_cpg_dict[contig_position] = read_pos

        # find location of ALL cpgs
        all_cpgs = re.finditer("[cC][gG]", temp_contig[1])
        for cpg in all_cpgs:
            # 2 Cs for each CpG (one on each strand)
            c_pos, g_pos = (cpg.start(0) + 1, cpg.start(0) + 2)
            read_pos = "NA"
            # check if position of C is in the dictionary of MspI CpGs
            if c_pos in MspI_cpg_dict:
                read_pos = MspI_cpg_dict[c_pos]
            outfile.write("%s\t%s\t%s\t%s\n" % (contig, c_pos,
                                                "+", read_pos))
            read_pos = "NA"
            if g_pos in MspI_cpg_dict:
                read_pos = MspI_cpg_dict[g_pos]
            outfile.write("%s\t%s\t%s\t%s\n" % (contig, g_pos,
                                                "-", read_pos))
        temp_contig = next(fasta)
    outfile.close()


@cluster_runnable
def mergeAndDrop(cpgs, infiles, outfile):
    '''merge bismark files with flat file of all cpgs
    into a single dataframe and drop unwanted columns'''
    # very memory intensive!

    # prepare dataframe with cpg locations
    all_cpgs = pd.io.parsers.read_csv(cpgs, sep='\t', comment='#')
    all_cpgs['location'] = (all_cpgs['contig'] + ":" +
                            all_cpgs['position'].apply(str))

    # for loop through files and merge each time with previous merged df
    for infile in infiles:
        sample_name = re.sub("_.*bismark.*", "", os.path.basename(infile))

        # prepare temporary dataframe
        column_names = ["contig", "start", "stop", "perc",
                        "meth", "unmeth"]
        temp_meth_df = pd.io.parsers.read_csv(infile, sep='\t',
                                              comment='#', header=None,
                                              names=column_names)
        temp_meth_df['location'] = (temp_meth_df['contig'] + ":" +
                                    temp_meth_df['start'].apply(str))
        temp_meth_df.drop(['contig', 'start', 'stop'], axis=1, inplace=True)
        temp_meth_df.columns = [hjoin([sample_name, "perc"]),
                                hjoin([sample_name, "meth"]),
                                hjoin([sample_name, "unmeth"]),
                                'location']

        # merge into temp dataframe
        merged = temp_meth_df.merge(all_cpgs, how="outer", left_on='location',
                                    right_on='location', sort=False)
        all_cpgs = merged.copy(deep=True)

    merged.drop('location', axis=1, inplace=True)
    merged.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def subsetToCovered(infile, outfile, cov_threshold=10):
    '''take a flat file containing methylation calls for all cpgs
    and output a flat file containing only sites above coverage threshold'''
    # to do: generalise function. Currently expects infile to contain meth
    # column followed by unmeth column for each sample

    raw = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    meth_cols = [x for x in raw.columns.tolist()
                 if re.match(".*meth|.*unmeth]", x)]
    meth_df = raw.ix[:, meth_cols]

    def minCov(row):
        # take every meth count column and add next column
        # (e.g unmeth count) and find min for row
        return min([np.nansum([i, k]) for i, k in zip(row[0::2], row[1::2])])

    meth_df['min_cov'] = meth_df.apply(minCov, axis=1)
    high_cov = raw[meth_df['min_cov'] >= cov_threshold]
    high_cov['cpgi'].fillna('Non-CpGIsland', inplace=True)
    high_cov.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def categorisePromoterCpGs(outfile, genome_fasta, annotations_database):
    '''extract promoter sequences and categorise them by CpG density'''

    def getGC(sequence):
        ''' calculate the O/E and GC content for a given sequence'''
        sequence = sequence.upper()
        C = sequence.count("C")
        G = sequence.count("G")
        seq_length = float(len(sequence)) - sequence.count("N")

        # sequence may be all Ns(!)
        if seq_length == 0 or C == 0 or G == 0:
            return None, None

        GC = (C + G) / seq_length

        OE = sequence.count("CG") * seq_length / float(C * G)

        return GC, OE

    def categorisePromoter(prom_seq, window=500, slide=1):
        ''' identify Low/Mid/High CpG promoters
        HCPs - At lease on 500bp interval with a GC fraction >=0.55 and O/E >= 0.6
        LCPs - No 500bp interval with CpG O/E >= 0.4
        MCPs - the rest'''
        HCP = False
        MCP = False
        LCP = False
        for seq_start in range(0, len(prom_seq) - window, slide):
            GC, OE = getGC(prom_seq[seq_start: seq_start + window])

            if GC >= 0.55 and OE >= 0.6:
                return "HCP"

            elif OE >= 0.4:
                MCP = True

            if GC:
                LCP = True

        if MCP:
            return "MCP"
        elif LCP:
            return "LCP"
        else:
            return None

    IxFA = IndexedFasta.PysamIndexedFasta(genome_fasta)
    IxFA.mConverter = nullConvert

    connect = sqlite3.connect(annotations_database)
    select_cmd = '''SELECT gene_id, contig, start, end, strand FROM
    geneset_all_gtf_genome_coordinates'''
    select = connect.execute(select_cmd)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("%s\n" % "\t".join(
            ("contig", "position", "feature", "CpG_density")))

        for promoter in select:
            gene_id, contig, start, end, strand = list(map(str, promoter))

            if strand == "+":
                if int(start) >= 2000:
                    prom_seq = IxFA.getSequence(
                        contig, strand, int(start) - 2000, int(start) + 500)
                    CpGs = getCpGs(prom_seq)
                    CpGs = [x + int(start) - 2000 for x in CpGs]

                else:
                    prom_seq = IxFA.getSequence(
                        contig, strand, 0, int(start) + 500)
                    CpGs = getCpGs(prom_seq)

            elif strand == "-":
                # if the gene ends within 2000 bp of the contig end,
                # this will skip the gene
                try:
                    prom_seq = IxFA.getSequence(
                        contig, strand, int(end) - 500, int(end) + 2000)
                except:
                    pass

                prom_seq = prom_seq[::-1]
                CpGs = getCpGs(prom_seq)

                CpGs = [int(end) + 2000 - x for x in CpGs[::-1]]

            prom_type = categorisePromoter(prom_seq)

            total_CG = len(CpGs)

            if total_CG == 0:
                CpG_density = "0"
            else:
                CpG_density = str(float(total_CG) / len(prom_seq))

            for CpG in CpGs:
                outf.write("%s\n" % "\t".join(
                    (contig, str(CpG), prom_type, CpG_density)))


@cluster_runnable
def findRepeatCpGs(outfile, genome_fasta, repeats_gff):
    '''extract repeats sequences and identify CpG locations'''

    IxFA = IndexedFasta.PysamIndexedFasta(genome_fasta)
    IxFA.mConverter = nullConvert

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("%s\n" % "\t".join(
            ("contig", "position", "feature", "CpG_density")))

        with IOTools.openFile(repeats_gff, "r") as inf:
            for line in inf:
                line = line.strip().split("\t")
                contig = line[0]
                start, stop = line[3:5]
                atr = line[8].split("; ")
                repeat_class = atr[0].split(" ")[1].replace('"', '')
                repeat_family = atr[1].split(" ")[1].replace('"', '')

                repeat_seq_f = IxFA.getSequence(
                    contig, "+", int(start), int(stop))

                CpGs_f = getCpGs(repeat_seq_f)

                CpGs = [x + int(start) for x in CpGs_f]
                total_CG = len(CpGs)

                CpGs.extend([x + 1 for x in CpGs])

                if total_CG == 0:
                    CpG_density = "0"
                else:
                    CpG_density = str(float(total_CG) / len(repeat_seq_f))

                for CpG in CpGs:
                    outf.write("%s\n" % "\t".join(
                        (contig, str(CpG), repeat_class, CpG_density)))

                    outf.write("%s\n" % "\t".join(
                        (contig, str(CpG), repeat_family, CpG_density)))


@cluster_runnable
def findCpGsFromBed(outfile, genome_fasta, bed, feature, both_strands=True):
    '''extract sequences for bed entries and identify CpG locations'''

    IxFA = IndexedFasta.PysamIndexedFasta(genome_fasta)
    IxFA.mConverter = nullConvert

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("%s\n" % "\t".join(
            ("contig", "position", "feature", "CpG_density")))

        with IOTools.openFile(bed, "r") as inf:
            for line in inf:
                line = line.strip().split("\t")
                contig, start, stop = line[0:3]

                seq = IxFA.getSequence(contig, "+", int(start), int(stop))

                CpGs = getCpGs(seq)

                CpGs = [x + int(start) for x in CpGs]

                seq_length = len(seq)

                if both_strands:
                    CpGs.extend([x + 1 for x in CpGs])
                    seq_length *= 2

                total_CG = len(CpGs)

                if total_CG == 0:
                    CpG_density = "0"
                else:
                    CpG_density = str(float(total_CG) / seq_length)

                for CpG in CpGs:
                    outf.write("%s\n" % "\t".join(
                        (contig, str(CpG), feature, CpG_density)))


@cluster_runnable
def mergeCpGAnnotations(meth_inf, prom_inf, repeat_inf, hcne_inf, dmr_inf,
                        outfile):
    '''merge together the CpG annotations for plotting'''

    df = pd.read_table(meth_inf, usecols=[
        "Germline-Dex-mean", "Germline-Veh-mean", "Liver-Dex-mean",
        "Liver-Veh-mean", "contig", "position"])
    df.columns = [x.replace("-mean", "") for x in df.columns]
    df['position'] = df['position'] - 1
    df.set_index(["contig", "position"], inplace=True)

    def readAndMerge(df, cpg_annotations_inf):
        tmp_df = pd.read_table(cpg_annotations_inf)
        tmp_df.set_index(["contig", "position"], inplace=True)
        tmp_df = tmp_df.join(df, how="inner")
        tmp_df = pd.melt(tmp_df, id_vars=["feature", "CpG_density"])
        return tmp_df

    df_promoter = readAndMerge(df, prom_inf)
    df_repeats = readAndMerge(df, repeat_inf)
    df_hcne = readAndMerge(df, hcne_inf)
    df_dmr = readAndMerge(df, dmr_inf)

    final_df = pd.concat([df_promoter, df_repeats, df_hcne, df_dmr])
    final_df.to_csv(outfile, sep="\t")


def plotCpGAnnotations(infile, outfile_hist, outfile_box):
    ''' make histogram and boxplots for the CpGs facetted per annotation'''
    df = pd.read_table(infile, sep="\t")
    r_df = pandas2ri.py2ri(df)

    plotter = r('''
    function(df){
    library(scales)
    library(plyr)
    library(ggplot2)

    # shouldn't hard code the features to retain
    retain_features = c("HCP", "MCP", "LCP", "HCNE", "DMR", "LTR",
                        "SINE", "Alu", "LINE", "ERVK")
    df = df[df$feature %%in%% retain_features,]
    df$feature = factor(df$feature, levels = retain_features)

    # count instances per facet
    df.n <- ddply(.data=df, .(variable, feature),
    summarize, n=paste("n =", length(value)))

    # we just want counts in the top row
    df.n = df.n[df.n$variable=="Germline-Dex",]

    ltxt = element_text(size = 20)
    mtxt = element_text(size = 15)
    mtxt90 = element_text(size = 15, angle=90, vjust=0.5, hjust=1)

    theme_custom = theme(aspect.ratio=1,
    axis.text.y = mtxt, axis.text.x = mtxt90,
    axis.title.y = ltxt, axis.title.x = ltxt,
    strip.text.x=ltxt, strip.text.y=mtxt,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(colour="grey60", fill="grey90"))

    p = ggplot(df, aes(value)) +
    geom_histogram(aes(y=100*(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
    binwidth=10, fill="midnightblue") +
    facet_grid(variable~feature) +
    geom_text(data=df.n, aes(x=60, y=90, label=n),
    colour="black", inherit.aes=FALSE, parse=FALSE, size=5) +
    theme_bw() +
    theme_custom +
    xlab("Methylation (%%)") +
    ylab("CpGs (%%)") +
    scale_x_continuous(breaks=c(0,50,100)) +
    scale_y_continuous(breaks=c(0,25, 50, 75, 100))

    ggsave("%(outfile_hist)s", width=16, height=8)

    # for the box plot we need to bin the density values
    df['binned_density'] <- .bincode(x = 100*df$CpG_density,
                                     breaks = c(0,1,3,5,7,10,1000))

    df['binned_density'] = revalue(as.factor(df$binned_density),
                                   c("1"="<1", "2"="1-3", "3"="3-5",
                                     "4"="5-7", "5"="7-10", "6"=">10"))

    p = ggplot(df, aes(as.factor(binned_density), value, group=binned_density)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(notch=TRUE, outlier.shape = NA, fill="grey80", width = 0.8) +
    stat_summary(fun.y="median", colour="red", geom="point", size=1) +
    coord_flip() +
    facet_grid(variable~feature) +
    theme_bw() +
    theme_custom +
    ylab("Methylation (%%)") +
    xlab("CpG Density (%%)") +
    scale_y_continuous(breaks=c(0,50,100))

    ggsave("%(outfile_box)s", width=16, height=8)

    }''' % locals())

    plotter(r_df)


def identifyExperimentalFactors(df):
    '''takes a data frame and returns the samples, treatment groups,
    replicate numbers and tissues based on tissue-treatment-replicate
    nomenclature'''

    columns = []
    for col in list(df.columns.values):
        m = re.match("\S*-\S*-\S*-perc", col)
        if m is not None:
            columns.append(col)

    samples = [re.sub("-perc", "", x) for x in columns]
    treatments = list(set([x.split("-")[1] for x in samples]))
    tissues = list(set([x.split("-")[0] for x in samples]))
    replicates = list(set([x.split("-")[2] for x in samples]))

    return samples, treatments, tissues, replicates


@cluster_runnable
def addTreatmentMean(infile, outfile):
    '''take a flat file containing fraction methylation per sample and
    calculate treatment means based on column names'''

    df = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    samples, treatments, tissues, replicates = identifyExperimentalFactors(df)

    for tissue in tissues:
        for treatment in treatments:
            p_columns = []
            for replicate in replicates:
                p_columns.append(hjoin([tissue, treatment, replicate, "perc"]))
            temp_df = df[p_columns]

            df[hjoin([tissue, treatment, "mean"])] = temp_df.apply(
                np.mean, axis=1)
            df[hjoin([tissue, treatment, "var"])] = temp_df.apply(
                np.var, axis=1)

    df.to_csv(outfile, sep="\t", index=False, na_rep="NA")


# this should be transferred into farm.py
@cluster_runnable
def splitDataframeClusters(infile, prefix, suffix):
    df = pd.read_csv(infile, comment='#', sep="\t")
    df.sort(columns=["contig", "position"], inplace=True)
    df.drop(["strand", "read_position"], inplace=True, axis=1)
    positions = df['position'].astype('int').tolist()
    contigs = df['contig'].tolist()
    # df.set_index('position', inplace=True)
    current_pos, count, splits, n, first_pos = (0,) * 5
    current_contig = ""
    size = 10
    distance = 100
    split_at = 10000
    for ix in range(0, len(positions)):
        next_pos = positions[ix]
        next_contig = contigs[ix]
        if (((next_pos < current_pos + distance) &
             (next_contig == current_contig))):
            count += 1
        else:
            if count >= size:
                temp_cluster_df = df[first_pos:ix + 1]
                if n >= split_at:
                    identifier = splits * split_at
                    filename = "%(prefix)s%(identifier)i%(suffix)s" % locals()
                    cluster_df.to_csv(
                        filename, index=False, header=True, sep="\t",
                        dtype={'position': int})
                    n = 0
                    splits += 1
                if n == 0:
                    cluster_df = temp_cluster_df
                    n = count
                else:
                    cluster_df = cluster_df.append(temp_cluster_df)
                    n += count
            current_pos = next_pos
            current_contig = next_contig
            first_pos = ix + 1
            count = 0

    # make sure final clusters are written out into smaller final file
    if n < split_at:
        identifier = splits * split_at
        filename = "%(prefix)s%(identifier)i%(suffix)s" % locals()
        cluster_df.to_csv(
            filename, index=False, header=True, sep="\t",
            dtype={'position': int})


@cluster_runnable
def calculateCoverage(infile, outfile):
    ''' calculate the coverage at CpG islands
    and non-CpG Islands '''

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("\t".join(("coverage", "cpgi")) + "\n")

        with IOTools.openFile(infile) as inf:

            header = inf.next().strip().split("\t")
            coverage_cols = ["meth" in x for x in header]
            n_samples = sum(coverage_cols) / 2
            cpgi_col = header.index('cpgi')

            line_n = 0
            for line in inf:
                if line_n % 1000000 == 0:
                    E.info("parsed %i lines" % line_n)

                line_n += 1
                line = line.strip().split("\t")

                cov_values = [line[x] for x in range(0, len(line))
                              if coverage_cols[x]]

                cov = np.nansum(np.genfromtxt(
                    np.array(cov_values))) / n_samples

                if line[cpgi_col] == "CpGIsland":
                    cpgi = "CpGIsland"
                else:
                    cpgi = "Non-CpGIsland"

                outf.write("\t".join(map(str, (cov, cpgi))) + "\n")


@cluster_runnable
def plotCoverage(infile, outfiles):
    ''' plot the coverage at CpG islands
    and non-CpG Islands '''

    bar_plot, pie_plot = outfiles
    binned_out = P.snip(bar_plot, "_bar.png") + "_binned.tsv"

    plotCoverage = r('''
    library(ggplot2)
    library(plyr)
    library(RColorBrewer)

    function(infile, binned_out, bar_plot, pie_plot){

    df = read.table(infile, sep="\t", header=TRUE)

    df$binned = .bincode(
        df$coverage, breaks = c(0,1,5,10,20,1000000), include.lowest=TRUE)

    df$binned = factor(df$binned)

    df$binned = revalue(df$binned,
        c("1"="< 1", "2"="[ >1 - 5 ]", "3"="[ >5 - 10 ]",
          "4"="[ >10 - 20 ]", "5"="> 20"))

    write.table(df, binned_out, quote=FALSE, sep="\t", row.names=FALSE)

    s_f = scale_fill_manual(name="Coverage", values=brewer.pal(5, "YlOrRd"))

    x = xlab("")
    y = ylab("Fraction")
    l_txt = element_text(size=10)
    m_txt = element_text(size=12)
    t = theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        aspect.ratio=1)

    t2 = theme(
        axis.text.x=m_txt,
        axis.text.y=l_txt,
        axis.title.y=l_txt,
        legend.title=l_txt,
        legend.text=l_txt,
        strip.text=l_txt,
        aspect.ratio=1)

    p = ggplot(df, aes(x=cpgi, fill=as.factor(binned))) +
        geom_bar(position="fill", stat="bin") +
        theme_bw() + s_f + x + t2 + y

    ggsave(bar_plot, height=4, width=5)

    p = ggplot(df, aes(factor(0), fill=as.factor(binned))) +
        geom_bar(position="fill", stat="bin", width=1) +
        coord_polar(theta = "y") +
        facet_wrap(~cpgi) +
        theme_bw() + s_f + x + t + y

    ggsave(pie_plot, height=3, width=6)
    }''')

    plotCoverage(infile, binned_out, bar_plot, pie_plot)


@cluster_runnable
def plotMethFrequency(infile, outfile):

    plotMethFrequency = r('''
    function(infile, outfile){
    library(ggplot2)
    library(grid)
    library(RColorBrewer)

    df = read.table(infile, sep="\t", header=TRUE)

    plotter = function(df, plotname, facet='both'){

        l_size = 30
        l_txt = element_text(size = l_size)

        p = ggplot(df, aes(x=treatment_replicate,
                            y=X0, fill=value)) +
        geom_bar(position="fill", stat="identity") +
        xlab("") +
        ylab("Fraction") +
        scale_fill_manual(name="Methylation (%)",
                          values = rev(colorRampPalette(
                                           brewer.pal(9, "RdBu"))(11))) +
        scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1)) +
        theme_bw() +
        theme(
        axis.text.x = element_text(angle=90, size=l_size, hjust=1, vjust=0.5),
        axis.text.y = l_txt,
        axis.title.y = l_txt,
        aspect.ratio=1,
        legend.text = l_txt,
        legend.title = l_txt,
        legend.key.size = unit(1.5, "cm"),
        strip.text = l_txt)

        if (facet=='both'){
        p = p + facet_grid(tissue~cpgi)
        ggsave(plotname, height=10, width=12)}

        else if (facet=='cpgi'){
        p = p + facet_grid(~cpgi)
        ggsave(plotname, height=7, width=12)}
        }

    plotter(df, outfile, facet='both')

    for(tissue in levels(factor(df$tissue))){
    tmp_df = df[df$tissue==tissue,]
    plotfile = gsub(".png", paste0("_", tissue, ".png"), outfile)
    plotter(tmp_df, plotfile, facet='cpgi')
    }}''')

    plotMethFrequency(infile, outfile)


@cluster_runnable
def calculateM3DStat(infile, outfile, design,
                     pair=None, groups=None):
    '''calculate M3D stats from dataframe'''
    # assumes meth columns will be named *-*-*-meth
    df = pd.read_csv(infile, sep="\t")
    design_df = pd.read_csv(design, sep="\t")

    design_df = design_df.ix[design_df['include'] == 1, ]

    if pair:
        design_df.set_index(design_df['track'], drop=False, inplace=True)
        pair1, pair2, = pair
        keep_tracks = [tracks for tracks in design_df['track'].tolist()
                       for track in (pair1, pair2) if track in tracks]
        design_df_subset = design_df.ix[keep_tracks, :]
        # create conditions from group and design_df
        conditions = [
            group for track in design_df_subset['track']
            for group in groups
            if group in track]
    else:
        # use the first pair specified in the design table
        pair = design_df['pair'].tolist()[0]
        design_df_subset = design_df[design_df['pair'] == pair]
        conditions = design_df_subset['group']

    samples = design_df_subset['track']
    keep_columns = ["contig", "position"]
    keep_columns.extend([x + "-perc" for x in samples])
    keep_columns.extend([x + y for x in design_df_subset['track']
                         for y in ["-unmeth", "-meth"]])

    df = df.ix[:, keep_columns]
    ncols = len(samples)
    samples = '","'.join(samples)
    conditions = '","'.join(conditions)
    r_df = pandas2ri.py2ri(df)

    out = open(outfile, "w")
    out.write("ncols: %s\n" % ncols)
    out.write("conditions: %s\n" % conditions)
    out.write("samples: %s\n" % samples)
    out.write("pair: %s\n" % pair)
    out.write("groups: %s\n" % groups)
    out.write("%i, %i" % df.shape)
    out.close()

    func = r('''
    function(df){
    library(M3D)
    library(BiSeq)
rowRanges <- GRanges(seqnames = Rle(as.character(df$contig)),
    ranges = IRanges(start = df$position, end= df$position))
    colData <- DataFrame(group = c("%(conditions)s"),
    row.names = c("%(samples)s"))
    methReads <- matrix(as.integer(unlist(df[,grep("*[.]meth", colnames(df),
    value=FALSE)])),ncol=%(ncols)i)
    unmethReads <- matrix(as.integer(unlist(df[,grep("*[.]unmeth",colnames(df),
    value=FALSE)])),ncol=%(ncols)i)
    totalReads <- methReads + unmethReads
    rawData=BSraw(rowRanges = rowRanges, colData = colData,
    totalReads = totalReads, methReads = methReads)
    clust.unlim <-clusterSites(object = rawData,  perc.samples = 1,
    min.sites = 10, max.dist = 100)
    clust.unlim_GR <- clusterSitesToGR(clust.unlim)
    overlaps<-findOverlaps(clust.unlim_GR,rowRanges(rawData))
    MMDlist<-M3D_Wrapper(rawData, overlaps)
    M3Dstat<-MMDlist$Full-MMDlist$Coverage
    adjusted_colnames <- gsub(" ", "_", colnames(M3Dstat))
    M3Dstat_df <- as.data.frame(M3Dstat)
    colnames(M3Dstat_df) <- adjusted_colnames
    ranges_df <- data.frame(seqnames=seqnames(clust.unlim_GR,
    start=start(clust.unlim_GR)-1,end=end(clust.unlim_GR))
    complete_df <- cbind(ranges_df, M3Dstat_df)
    write.table(complete_df,file="%(outfile)s",sep="\t",quote=F,row.names=F)
    return(M3Dstat)
    }''' % locals())

    df.to_csv(outfile + ".pd", sep="\t")
    # out.close()
    M3Dstat = func(r_df)


def M3Dstat2pvalue(df, columns, pair):

    stats = importr('stats')

    def meltAndPivot(df, columns):
        melt = pd.melt(df, id_vars=columns)
        melt['value'] = melt['value'].astype(float)
        pivot = pd.DataFrame(melt.pivot_table(
            index=columns, values="value", aggfunc=np.mean))
        return pivot

    def sampler(M3Dstat, vector, n):
        rand = np.random.choice(vector, size=n, replace=True, )
        return len(rand[rand > M3Dstat]) / n

    def addPvalues(df, within_df, repeat):
        df['p_value'] = df['value'].apply(
            sampler, args=(within_df['value'], repeat))
        df['p_value_adj'] = stats.p_adjust(FloatVector(
            df['p_value']), method='BH')

    E.info("dimensions:")
    E.info(df.shape)

    M3D_within_col = []
    M3D_between_col = []

    # create regexes from pair in both directions to find columns
    for p in [pair, pair[::-1]]:
        regex = "%s-\d__vs__%s-\d" % (p[0], p[0])
        M3D_within_col.extend([col for col in df.columns
                               if re.search(regex, col)])
        regex = "%s-\d__vs__%s-\d" % (p[0], p[1])
        M3D_between_col.extend([col for col in df.columns
                                if re.search(regex, col)])

    within = copy.copy(columns)
    within.extend(M3D_within_col)
    between = copy.copy(columns)
    between.extend(M3D_between_col)

    between_pivot = meltAndPivot(df.ix[:, between], columns)
    within_pivot = meltAndPivot(df.ix[:, within], columns)

    addPvalues(between_pivot, within_pivot, 100000.0)
    addPvalues(within_pivot, within_pivot, 100000.0)

    between_pivot.reset_index(inplace=True)
    within_pivot.reset_index(inplace=True)

    return between_pivot, within_pivot


@cluster_runnable
def calculateM3DSpikepvalue(infiles, outfile, design):
    '''takes a list of M3D stats files and a design dataframe and outputs
    p-values for each cluster for the pairwise comparison'''

    outfile_within = P.snip(outfile, "_between.tsv") + "_within.tsv"

    design_df = pd.read_csv(design, sep="\t")
    design_df = design_df.ix[design_df['include'] == 1, ]

    first_pair = design_df['pair'].tolist()[0]
    design_df_subset = design_df[design_df['pair'] == first_pair]

    # need to create list along the lines of:
    # ["Germline-Saline", "Germline-Dex"]
    # here, this is acheived by removing replicates from track names and
    # outputting unique values - should use PipelineTracks.py
    pair = list(set([re.sub("-(\d+)", "", x)
                     for x in design_df_subset['track'].tolist()]))

    df = pd.read_csv(infiles[0], sep="\t")
    if len(infiles) > 0:
        for next_file in infiles[1:]:
            temp_df = pd.read_csv(next_file, sep="\t")
            df = df.append(temp_df)
    df.set_index(df['seqnames'], inplace=True)

    cluster_characteristics = ["contig", "initial", "change",
                               "size", "index", "id"]
    seqnames_split = df['seqnames'].apply(lambda x: x.split("_"))
    seqnames_split = [x for x in seqnames_split]

    cluster_values_df = pd.DataFrame(data=seqnames_split,
                                     index=df.index)
    cluster_values_df.columns = cluster_characteristics
    concat_df = pd.concat([df, cluster_values_df], axis=1)
    # remove after testing
    concat_df.to_csv(outfile, index=False, header=True, sep="\t")

    between, within = M3Dstat2pvalue(concat_df, cluster_characteristics,
                                     pair)

    between.to_csv(outfile, index=False, header=True, sep="\t")
    within.to_csv(outfile_within, index=False, header=True, sep="\t")

    r_between_melt = pandas2ri.py2ri(between)
    r_within_melt = pandas2ri.py2ri(within)

    # move plotting to reporting. output an aggregated df.
    base = "power.dir/M3D"
    plotter = r('''
    library(ggplot2)
    function(df){
    agg = aggregate(df$p_value_adj, by=list(df$size, df$change),
    FUN=function(x) length(x[x<0.05])/length(x))
    colnames(agg) = c("size","change","power")

    text_theme = element_text(size=20)
    text_theme_s = element_text(size=15)
    text_theme_xs = element_text(size=10)
    text_theme_s_x = element_text(size=15, angle=90)
    axis_title_x = element_text(size=15, vjust=-2)
    axis_title_y = element_text(size=15, vjust=2)
    text_theme_l = element_text(size=25)

    plot_theme = theme(
    axis.title.x = axis_title_x,
    axis.title.y = axis_title_y,
    axis.text.x = text_theme_s_x,
    axis.text.y = text_theme,
    legend.text = text_theme_s,
    legend.title = text_theme_s,
    title = text_theme_s,
    aspect.ratio=1,
    panel.grid.minor = element_blank())

    p = ggplot(agg,aes(x=as.numeric(change), y=as.numeric(as.character(size)),
    fill=as.numeric(power))) +
    geom_tile() +
    scale_fill_continuous(limits=c(0,1), name = "Power",
                          low="grey95", high="mediumpurple4") +
    scale_y_continuous(breaks=seq(1,9,1)) +
    scale_x_continuous(breaks=seq(-50,50,10), limits=c(-50,50)) +
    theme_bw() +
    plot_theme +
    ylab("Size of spike-in region") +
    xlab("Methylation change (%%)")
    ggsave(paste0('%(base)s',"_change_vs_size_power_heatmap.png"),
    plot = p, width = 5, height = 5)

    p = ggplot(agg,aes(x=as.numeric(change),
    y=as.numeric(power),col=factor(as.numeric(as.character(size))))) +
    geom_line(size=2) + scale_colour_discrete(name = "Size") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    scale_x_continuous(breaks=seq(-50,50,10), limits=c(-50,50)) +
    theme_bw() +
    plot_theme +
    xlab("Methylation change (%%)") +
    ylab("Power")
    ggsave(paste0('%(base)s',"_change_vs_size_power_lineplot.png"),
    plot = p, width = 5, height = 5)

    if (sum(agg$power>0.8 & agg$change>0) > 0 &
        sum(agg$power>0.8 & agg$change<0) >0){
      powered_df = (agg[agg$power>0.8 & agg$change>0,])
      powered_df2 = (agg[agg$power>0.8 & agg$change<0,])

      powered_agg = aggregate(as.numeric(powered_df$change),
          by=list(powered_df$size),FUN=min)
      powered_agg2 = aggregate(as.numeric(powered_df2$change),
          by=list(powered_df2$size),FUN=max)
      powered_agg['direction'] = "gain"
      powered_agg2['direction'] = "loss"
      merge = rbind(powered_agg, powered_agg2)
      colnames(merge) = c("size","change","direction")
      merge$size = as.numeric(as.character(merge$size))
      merge$change = as.numeric(as.character(merge$change))

      p = ggplot(merge, aes(x=size, y=change, group=direction)) +
      geom_line(size=3) +
      theme_bw() + plot_theme +
      scale_y_continuous(breaks=seq(-50,50,10), limits=c(-50,50)) +
      scale_x_continuous(breaks=seq(1,9,1), limits=c(1,8)) +
      xlab("Size") + ylab("Methylation difference") +

      ggtitle ("M3D - Mimimum methylation\n difference to reach 0.8 power")
      ggsave(paste0('%(base)s',"_powered.png"),
      plot = p, width = 5, height = 5)}}
    ''' % locals())

    plotter(r_between_melt)


@cluster_runnable
def calculateM3Dpvalue(infiles, outfile, pair):
    '''takes a list of M3D stats files and a design dataframe and outputs
    p-values for each cluster for the pairwise comparison'''
    # to do: use design dataframe to identify groups
    # currently hardcoded in pipeline

    outfile_within = P.snip(outfile, "_between.tsv") + "_within.tsv"
    plots_dir = os.path.dirname(outfile)

    df = pd.read_csv(infiles[0], sep="\t")
    if len(infiles) > 0:
        for next_file in infiles[1:]:
            temp_df = pd.read_csv(next_file, sep="\t")
            df = df.append(temp_df)
    df.set_index(df['seqnames'], inplace=True)
    between, within = M3Dstat2pvalue(
        df, ["seqnames", "start", "end"], pair)

    between.to_csv(outfile, index=False, header=True, sep="\t")
    within.to_csv(outfile_within, index=False, header=True, sep="\t")

    r_between_melt = pandas2ri.py2ri(between)
    r_within_melt = pandas2ri.py2ri(within)

    # move plotting to report. i.e. output a merged dataframe with a
    # "cluster" column
    title = "_vs_".join(pair)
    base = "%(plots_dir)s/M3D_%(title)s" % locals()
    plotter = r('''
    library(ggplot2)
    function(df, df2){
    df["cluster"] = "between"
    df2["cluster"] = "within"
    merge = rbind(df,df2)
    base_theme = theme(axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15), legend.text=element_text(size=15),
    axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
    legend.title=element_text(size=15), legend.position=c(0.9,0.9))

    p = ggplot(merge, aes(x=p_value)) + geom_density(aes(colour=cluster)) +
        base_theme + scale_colour_discrete(name="Comparison") +
        xlab("p-value (unadjusted)") + ggtitle('%(title)s')
    p2 = ggplot(merge, aes(x=value)) + geom_histogram(aes(fill=cluster),
         binwidth=0.05, position="identity",alpha=0.5) + base_theme +
         scale_fill_discrete(name="Comparison") + xlab("M3D statistic") +
         ggtitle('%(title)s')
    p3 = ggplot(merge, aes(x=value)) +
         stat_density(aes(fill=cluster),
                      position="identity",alpha=0.5) + base_theme +
         xlab("M3D statistic") + scale_fill_discrete(name="Comparison")+
         ggtitle('%(title)s')
    p4 = ggplot(merge, aes(x=value)) + stat_density() +
         xlab("M3D statistic") + facet_grid(cluster~.) + base_theme +
         ggtitle('%(title)s')

    ggsave(paste0('%(base)s', "_p-value_dist.png"),
    plot = p, width = 5, height = 5)
    ggsave(paste0('%(base)s', "_stat_hist.png"),
    plot = p2, width = 5, height = 5)
    ggsave(paste0('%(base)s', "_stat_dist.png"),
    plot = p3, width = 5, height = 5)
    ggsave(paste0('%(base)s', "_stat_dist_facet.png"),
    plot = p4, width = 5, height = 5)}''' % locals())

    plotter(r_between_melt, r_within_melt)


@cluster_runnable
def summariseM3D(infile, outfile, threshold):
    '''count number of clusters passing threshold'''
    out = open(outfile, "w")
    df = pd.read_csv(infile, sep="\t")
    out.write("group1\tgroup2\ttotal\tsignificant\n")
    # remove any rows without 'value' column entry
    df_complete = df.ix[df['value'] > -1, ]
    total = len(df_complete.iloc[:, 1])
    g1, g2 = os.path.basename(infile).split("_vs_")
    g2 = g2.split("_")[0]
    df_thresh = df_complete.ix[df_complete['p_value_adj'] < threshold, ]
    significant = len(df_thresh.iloc[:, 1])
    out.write("%(g1)s\t%(g2)s\t%(total)s\t%(significant)s\n" % locals())
    out.close()


@cluster_runnable
def calculateBiSeqStat(infile, outfile):
    '''calculate BiSeq stats from dataframe'''
    # assumes meth columns will be named *-*-*-meth
    df = pd.read_csv(infile, sep="\t")
    r_df = pandas2ri.py2ri(df)
    samples = [re.sub(".perc", "", x) for x in df.columns
               if re.match(".*perc", x)]
    conditions = [x.split("-")[1] for x in samples]
    samples = '","'.join(samples)
    conditions = '","'.join(conditions)
    func = r('''library(BiSeq)
    library(M3D)
    function(df){
    rowData <- GRanges(seqnames = Rle(as.character(df$contig)),
    ranges = IRanges(start = df$position, end= df$position))
    colData <- DataFrame(group = c("%(conditions)s"),
    row.names = c("%(samples)s"))
    methReads <- matrix(as.integer(unlist(df[,grep("*[.]meth", colnames(df),
    value=FALSE)])),ncol=6)
    unmethReads <- matrix(as.integer(unlist(df[,grep("*[.]unmeth",colnames(df),
    value=FALSE)])),ncol=6)
    totalReads <- methReads + unmethReads
    rawData=BSraw(rowData = rowData, colData = colData,
    totalReads = totalReads, methReads = methReads)
    relData <- rawToRel(rawData)
    print(rawData)
    clust.unlim <-clusterSites(object = rawData,  perc.samples = 1,
    min.sites = 10, max.dist = 100, mc.cores=6)
    ind.cov <- totalReads(clust.unlim) > 0
    quant <- quantile(totalReads(clust.unlim)[ind.cov], 0.95)
    clust.lim <- limitCov(clust.unlim, maxCov = quant)
    predictedMeth <- predictMeth(object = clust.lim, mc.cores=6)
    betaResults <- betaRegression(formula = ~group,link = "probit",
    object = predictedMeth,type = "BR",mc.cores=6)
    predictedMethNull <- predictedMeth
    group_levels <- levels(colData(predictedMethNull)$group)
    null_array=NULL
    last_g1=0
    last_g2=1
    n=1
    for (x in colData(predictedMethNull)$group){
    if (x == group_levels[1])
    {null_array[n] = last_g1
    last_g1 = abs(last_g1-1)
    n = n+1}
    else
    {null_array[n] = last_g2
    last_g2 = abs(last_g2-1)
    n = n+1}}
    colData(predictedMethNull)$group.null <- null_array
    betaResultsNull <- betaRegression(formula = ~group.null,
    link = "probit",object = predictedMethNull, type="BR",
    mc.cores=6)
    vario <- makeVariogram(betaResultsNull)
    auto_sill = mean(vario$variogram$v[100:
    length(vario$variogram$v[,2]),2])
    vario.sm <- smoothVariogram(vario, sill = auto_sill)
    vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
    vario.sm$pValsList <- vario.aux$pValsList
    locCor <- estLocCor(vario.sm)

    clust.unlim_GR <- clusterSitesToGR(clust.unlim)
    ranges_df <- data.frame(seqnames=seqnames(clust.unlim_GR),
    start=start(clust.unlim_GR)-1,end=end(clust.unlim_GR))
    complete_df <- cbind(ranges_df, M3Dstat_df)
    write.table(complete_df,file="%(outfile)s",sep="\t",quote=F)
    return(M3Dstat)
    }''' % locals())
    BiSeqstat = func(r_df)


@cluster_runnable
def summaryPlots(infile, outfile):
    '''take a flat file containing average methylation values and
    produce various summary plots

    Plots:
    1. Read position vs. methylation level
    2. Distributions of average methylation, split by group and CpGI/Non-CpGI
    3. Methylation distributions by sample, split by CpGI/Non-CpGI
    4. Heatmap of Spearman's rho split by CpGI/Non-CpGI
    5. Correlation between mean value pairs for each treatment pair,
       split by CpGI/Non-CpGI
    6. Read position vs. methylation difference for each treatment pair

    to do:
    1. Add all vs. all correlations plots
    test!
    '''

    df = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    samples, treatments, tissues, replicates = identifyExperimentalFactors(df)

    header = [re.sub("-", "_", x) for x in df.columns.values.tolist()]
    df.columns = header

    r_dataframe = pandas2ri.py2ri(df)

    base = outfile
    gr = importr('grDevices')
    out = open(outfile + ".log", "w")

    plotFuncReadPositionMethylation = r("""
    library(ggplot2)
    library(reshape)
    function(df){
    read_length = max(df$read_position[is.finite(df$read_position)])
    df2 = df[, grepl("mean", names(df))]
    df_agg = aggregate(df2, by=list(df$read_position),FUN=mean, na.rm=TRUE)
    df_agg_m = melt(df_agg,id="Group.1")
    p = ggplot(df_agg_m, aes(x=Group.1,y = value)) +
    geom_line(aes(col=variable)) +
    geom_point(aes(col=variable), size=1.5) +
    ggtitle("Read Position vs. Methylation") +
    scale_colour_discrete(name = "Sample") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5))+
    xlab("Read position") +
    ylab("Methylation by read position")

    ggsave(paste0('%(base)s','-read_position_meth.png'),
    plot = p, width = 5, height = 5)
    }""" % locals())

    plotFuncReadPositionMethylation(r_dataframe)
    out.write("plot saved at: %(base)s-read_position_meth.png\n" % locals())

    plotMethBySampleMean = r("""
    library(ggplot2)
    library(reshape)
    function(df){
    df2 = df[df$cpgi != "CpGIsland",]
    df3 = df2[, grepl("mean", names(df))]
    df4 = melt(df3)

    p=ggplot(df4,
    aes(y=as.numeric(value), x=variable)) +
    geom_violin(aes(col=variable)) +
    scale_colour_discrete(name = "Sample", guide = FALSE) +
    ylab("Fraction Methylated") + xlab("") +
    scale_x_discrete(labels = gsub("_", "-", levels(df4$variable))) +
    theme(legend.position="bottom",
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5)) +
    ggtitle("Methylation by treatment group: Non-CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_type_non_cpgi.png'),
    plot = p, width = 5, height = 5)

    df5 = df[df$cpgi == "CpGIsland",]
    df6 = df5[, grepl("mean", names(df))]
    df7 = melt(df6)

    p2=ggplot(df7,
    aes(y=as.numeric(value), x=variable)) +
    geom_violin(aes(col=variable)) +
    scale_colour_discrete(name = "Sample", guide = FALSE) +
    ylab("Fraction Methylated") + xlab("") +
    theme(legend.position="bottom",
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5)) +
    ggtitle("Methylation by treatment group: CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_type_cpgi.png'),
    plot = p2, width = 5, height = 5)
    }""" % locals())

    plotMethBySampleMean(r_dataframe)
    out.write("plot saved at: \
    %(base)s-methylation_by_sample_type_non_cpgi.png\n" % locals())
    out.write("plot saved at: \
    %(base)s-methylation_by_sample_type_cpgi.png\n" % locals())

    meth_freq_columns = ["_".join((x, y, z, "perc"))
                         for x in tissues
                         for y in treatments
                         for z in replicates]
    meth_freq_columns.append("cpgi")

    df_meth_freq = df[meth_freq_columns]
    header = [re.sub("_perc", "", x) for x in df.columns.values.tolist()]
    df.columns = header

    r_dataframe_meth_freq = pandas2ri.py2ri(df_meth_freq)

    # re-write with subfunctions in r function to reduce length of code
    # each plot is repeated twice!
    drop_cpgi = '!(colnames(df) %in% c("cpgi"))'
    plotMethBySampleAndHeatmap = r("""
    library(ggplot2)
    library(reshape)
    library(amap)
    function(df){
    df2 = df[, %(drop_cpgi)s]
    var = apply(df2,1,var)
    var[is.na(var)] <- 0
    mean = apply(df2,1,mean)
    df3 = df[(var > 0 & mean >= 0.1 & mean <= 0.9),]
    df4 = df3[df3$cpgi != "CpGIsland", %(drop_cpgi)s]
    df4_m = melt(df4)

    samples=data.frame(matrix(unlist(strsplit(as.character(df4_m$variable),"_")),
    ncol=4,byrow=TRUE))
    colnames(samples) = c("tissue","treatment","replicate","suffix")
    samples = samples[,1:3]
    df4_m = cbind(samples, df4_m)

    p=ggplot(df4_m,
    aes(x=as.numeric(value), group=variable)) +
    geom_density(aes(col=paste0(tissue,"-",treatment), linetype=replicate)) +
    scale_colour_discrete(name = "Treatment group") +
    scale_linetype_discrete(name = "Replicate") +
    xlab("Fraction Methylation") +
    ggtitle("Methylation by sample: Non-CpG Islands")

    #ggsave(paste0('%(base)s','_methylation_by_sample_non_cpgi.png'),
    #plot = p, width = 5, height = 5)

    df5 = as.matrix(sapply(df4[,apply(df4,2,max)>0], as.numeric))
    df6 = cor(df4, method="spearman")
    df6_m = melt(df6)
    p = ggplot(df6_m, aes(x=X1, y=X2, fill=value))+
    geom_tile() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_continuous(name = "rho") +
    ggtitle("non-CpG Islands")

    ggsave(paste0('%(base)s','_spearmans_heatmap_non_cpgi.png'),
    plot = p, width = 5, height = 5)

    df7 = df3[df3$cpgi == "CpGIsland", %(drop_cpgi)s]
    df7_m = melt(df7)
    samples=data.frame(matrix(unlist(strsplit(as.character(df7_m$variable),"_")),
    ncol=4,byrow=TRUE))
    colnames(samples) = c("tissue","treatment","replicate","suffix")
    samples = samples[,1:3]
    df7_m = cbind(samples, df7_m)

    p = ggplot(df7_m,
    aes(x=as.numeric(value), group=variable)) +
    geom_density(aes(col=paste0(tissue,"-",treatment), linetype=replicate)) +
    scale_colour_discrete(name = "Treatment group") +
    scale_linetype_discrete(name = "Replicate") +
    xlab("Fraction Methylation") +
    ggtitle("Methylation by sample: CpG Islands")

    #ggsave(paste0('%(base)s','_methylation_by_sample_cpgi.png'),
    #plot = p, width = 5, height = 5)

    df8 = as.matrix(sapply(df7[,apply(df7,2,max)>0], as.numeric))
    df9 = cor(df8, method="spearman")
    df9_m = melt(df9)
    write.table(file="/ifs/projects/proj034/test.tsv", df7, sep="\t")
    write.table(file="/ifs/projects/proj034/test2.tsv", df8, sep="\t")
    p2 = ggplot(df9_m, aes(x=X1, y=X2, fill=value))+
    theme(axis.text.x = element_text(angle = 90)) +
    geom_tile()  + xlab("") + ylab("") +
    scale_fill_continuous(name = "rho")  +
    ggtitle("CpG Islands")

    ggsave(paste0('%(base)s','_spearmans_heatmap_cpgi.png'),
    plot = p2, width = 5, height = 5)

    library(pvclust)
    # redefine distance function to use spearman's rank correlation
    dist.pvclust <- function(x){
    res <- as.dist(1 - cor(x, method="spearman", use="pairwise.complete.obs"))
    attr(res,"method") <- "correlation"
    return(res)
    x <- as.matrix(x)
    P <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)}

    result <- pvclust(as.matrix(df5), method.hclust="ward", nboot=100)
    png(paste0('%(base)s','_hierarchical_clustering_non_cpgi.png'))
    par(mar=c(0, 4, 4, 2))
    plot(result, xlab="", sub="", print.num=F, cex.pv=0.7,
    col.pv=c(0,3,0), main="Non-CpG islands")
    dev.off()

    result2 <- pvclust(as.matrix(df8), method.hclust="ward", nboot=100)
    png(paste0('%(base)s','_hierarchical_clustering_cpgi.png'))
    par(mar=c(0, 4, 4, 2))
    plot(result2, xlab="", sub="", print.num=F, cex.pv=0.7,
    col.pv=c(0,3,0), main="CpG islands")
    dev.off()

    }""" % locals())

    plotMethBySampleAndHeatmap(r_dataframe_meth_freq)
    out.write("plot saved at: \
    %(base)s-_spearmans_heatmap_non_cpgi.png\n" % locals())
    out.write("plot saved at: \
    %(base)s-_spearmans_heatmap_cpgi.png\n" % locals())

    for tissue in tissues:
        for treatment_pair in itertools.combinations(treatments, 2):
            treatment1, treatment2 = list(map(str, treatment_pair))
            tissue = str(tissue)
            out.write("plotting tissue: %s, treatment: %s vs. %s\n"
                      % (tissue, treatment1, treatment2))

            label_base = "-".join([tissue, treatment1, treatment2])

            title = "Global Correlation: \n\
            %(tissue)s: %(treatment1)s vs. %(treatment2)s" % locals()

            # re-write this function to include r plotting function and then
            # pass two dfs in turn split by cpgi so that only one python
            # function is required
            plotFuncGlobalCorrelation = r("""
            library(ggplot2)
            function(df){
            df_r1rm = df[df$read_position!=1 | !(is.finite(df$read_position)),]
            p = ggplot(df_r1rm, aes(x=%(tissue)s_%(treatment1)s_mean,
            y=%(tissue)s_%(treatment2)s_mean)) +
            geom_point(size=0.5, alpha = 0.2) +
            xlab("%(treatment1)s") + ylab("%(treatment2)s") +
            theme(axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5)) +
            ggtitle('%(title)s') + facet_grid(.~cpgi)
            ggsave(paste0('%(base)s','-global_correlation',
            '%(label_base)s','.png'), plot = p, width = 5, height = 5)
            }""" % locals())

            plotFuncGlobalCorrelation(r_dataframe)
            out.write("plot saved at: \
            %(base)s-global_correlation_cpgi_%(label_base)s.png\n" % locals())

            title = "Read Position vs. Difference in Methylation:\n\
            %(tissue)s: %(treatment1)s vs. %(treatment2)s" % locals()
            plotFuncReadPositionDifference = r("""
            library(ggplot2)
            function(df){

            f <- function(x) {
            r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
            names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
            r}
            quantiles = quantile(df$%(tissue)s_%(treatment1)s_mean-
            df$%(tissue)s_%(treatment2)s_mean, probs = c(0.01,0.99))

            df$read_position = factor(df$read_position)
            p = ggplot(df, aes(x=factor(read_position),
            y = (%(tissue)s_%(treatment1)s_mean-\
            %(tissue)s_%(treatment2)s_mean))) +
            stat_summary(fun.data=f, geom="boxplot") + ggtitle('%(title)s') +
            theme(axis.text.x = element_text(size=5, angle=90)) +
            xlab("Read position") + ylim(quantiles) +
            ylab("Difference in methylation:\
            %(treatment1)s vs. %(treatment2)s")
            ggsave(paste0('%(base)s','-read_position_meth_diff_',
            '%(label_base)s','.png'), plot = p, width = 5, height = 5)
            }""" % locals())
            # note:need to check whether final plot function works:pep8 changes
            plotFuncReadPositionDifference(r_dataframe)
            out.write("plot saved at: \
            %(base)s-read_position_meth_diff_%(label_base)s.png\n" % locals())

    out.close()


@cluster_runnable
def spikeInClustersAnalysis(infile, outfile):
    infiles = open(infile, "r").readlines()
    samples = [re.sub("_10_pipeline.*", "", os.path.basename(x).strip())
               for x in infiles]
    samples_str = '","'.join(samples)
    outfile_log = open(outfile + ".log", "w")
    outfile_log.close()
    # outfile_log.write("samples: %s\n" % samples)
    samples_treatment = [x[1] for x in samples]
    samples_treatment = '","'.join(samples_treatment)
    # outfile_log.write("sample treatments: %s\n" % samples_treatment)
    infiles = '","'.join([x.strip() for x in infiles])
    # outfile_log.write("infiles: %s \n" % infiles)
    M3D = r('''library(BiSeq)
            library(M3D)
            files <- c("%(infiles)s")
            samples <- data.frame(row.names=c("%(samples_str)s"),
                       group=c("%(samples_treatment)s"))
            function(){
            rawData <-readBismark(files,samples)
            clust.unlim <- clusterSites(object = rawData, perc.samples = 6/6,
                                        min.sites = 15, max.dist = 100)
            clust.unlim_GR=clusterSitesToGR(clust.unlim)
            overlaps<-findOverlaps(clust.unlim_GR,rowData(rawData))
            MMDlist<-M3D_Wrapper(rawData, overlaps)
            M3Dstat<-MMDlist$Full-MMDlist$Coverage
            ranges_df <- data.frame(seqnames=seqnames(clust.unlim_GR),
                start=start(clust.unlim_GR)-1,end=end(clust.unlim_GR))
            M3Dstat_df <- as.data.frame(M3Dstat)
            adjusted_colnames = gsub(" ", "_", colnames(M3Dstat_df))
            colnames(M3Dstat_df) = adjusted_colnames
            complete_df = cbind(ranges_df, M3Dstat_df)
            write.table(complete_df,file="%(outfile)s",sep="\t",quote=F)}'''
            % locals())
    M3D()  # run r code
    # outfile_log.close()


@cluster_runnable
def spikeInClustersAnalysisBiSeq(infile, outfile):
    infiles = open(infile, "r").readlines()
    samples = [re.sub("_10_pipeline.*", "", os.path.basename(x).strip())
               for x in infiles]
    samples_str = '","'.join(samples)
    print("samples: %s\n" % samples)
    samples_treatment = [x[1] for x in samples]
    group1, group2 = (samples_treatment[0], samples_treatment[1])
    samples_treatment = '","'.join(samples_treatment)
    print("sample treatments: %s\n" % samples_treatment)

    infiles = '","'.join([x.strip() for x in infiles])
    print(infiles)
    base = P.snip(outfile, ".out")

    BiSeq_power_analysis = r('''
            library(ggplot2)
            sink("%(outfile)s")
            function(){
            text_theme = element_text(size=20)
            text_theme_s = element_text(size=15)
            text_theme_l = element_text(size=25)
            plot_theme = theme(axis.title.x = text_theme,
            axis.title.y = text_theme, axis.text.x = text_theme_s,
            axis.text.y = text_theme_s, legend.text = text_theme_s,
            legend.title = text_theme_s, title = text_theme_l)

            library(BiSeq)
            files <- c("%(infiles)s")
            samples <- data.frame(row.names=c("%(samples_str)s"),
                       group=c("%(samples_treatment)s"))

            rawData <-readBismark(files,samples)
            clust.unlim <- clusterSites(object = rawData, perc.samples = 6/6,
                                        min.sites = 15, max.dist = 100)
            ind.cov <- totalReads(clust.unlim) > 0
            quant <- quantile(totalReads(clust.unlim)[ind.cov], 0.95)
            clust.lim <- limitCov(clust.unlim, maxCov = quant)
            predictedMeth <- predictMeth(object = clust.lim, mc.cores=6)
            Group1 <- predictedMeth[,
                      colData(predictedMeth)$group == "%(group1)s"]
            Group2 <- predictedMeth[,
                      colData(predictedMeth)$group == "%(group2)s"]
            #mean.group1 <- rowMeans(methLevel(Group1))
            #mean.group2 <- rowMeans(methLevel(Group2))
            #mean_df = data.frame(mean.group1,mean.group2)
            #p = ggplot(mean_df,aes(x=mean.group1, y=mean.group2)) +
            #geom_point(alpha=0.25,size=1.5) +
            #plot_theme + xlab("Mean %(group1)s") + ylab("Mean %(group2)s")
            #ggsave(paste0('%(base)s',"_scatter_means.png"),
            #plot = p, width = 5, height = 5)

            betaResults <- betaRegression(formula = ~group,link = "probit",
                object = predictedMeth,type = "BR",mc.cores=6)
            predictedMethNull <- predictedMeth
            group_levels <- levels(colData(predictedMethNull)$group)
            null_array=NULL
            last_g1=0
            last_g2=1
            n=1
            for (x in colData(predictedMethNull)$group){
                if (x == group_levels[1])
                    {null_array[n] = last_g1
                    last_g1 = abs(last_g1-1)
                    n = n+1}
                else
                    {null_array[n] = last_g2
                    last_g2 = abs(last_g2-1)
                    n = n+1}}
            colData(predictedMethNull)$group.null <- null_array
            betaResultsNull <- betaRegression(formula = ~group.null,
                link = "probit",object = predictedMethNull, type="BR",
                mc.cores=6)
            vario <- makeVariogram(betaResultsNull)
            auto_sill = mean(vario$variogram$v[100:
                             length(vario$variogram$v[,2]),2])
            vario.sm <- smoothVariogram(vario, sill = auto_sill)
            vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
            vario.sm$pValsList <- vario.aux$pValsList
            locCor <- estLocCor(vario.sm)
            clusters.rej <- testClusters(locCor,FDR.cluster = 0.1)
            clusters.rej$clusters.reject
            clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = 0.05)
            DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)

            all_df = as.data.frame(do.call(rbind,
                (strsplit(levels(seqnames(rowData(rawData))),"_"))))
            dmr_df = as.data.frame(do.call(rbind,
                (strsplit(levels(factor(as.data.frame(seqnames(DMRs))[,1])),"_"))))
            col_n = c("contig","size","start","end","initial","change")
            colnames(dmr_df) = col_n
            colnames(all_df) = col_n
            dmr_df$DMR=1
            all_df$DMR=0
            cat_df = rbind(all_df,dmr_df)
            cat_df$change = gsub("m", "-", cat_df$change)
            agg_DMR= aggregate(cat_df$DMR,by=list(cat_df$size, cat_df$change),
                FUN=function(x) length(x[x==1])/length(x))
            colnames(agg_DMR) = c("size","change","power")
            p = ggplot(agg_DMR,aes(y=as.numeric(change),
            x=as.numeric(as.character(size)),fill=as.numeric(power))) +
            geom_tile() +
            scale_fill_continuous(limits=c(0,1), name = "Power") +
            scale_x_continuous(breaks=seq(1,11,2)) +
            scale_y_continuous(breaks=seq(-10,50,10), limits=c(-10,50)) +
            plot_theme + xlab("Size") + ylab("Change") + ggtitle ("BiSeq")
            ggsave(paste0('%(base)s',"_change_vs_size_power_heatmap.png"),
            plot = p, width = 5, height = 5)

            p = ggplot(agg_DMR,aes(x=as.numeric(change),
            y=as.numeric(power),col=factor(as.character(size)))) +
            geom_line(size=2) +
            scale_colour_discrete(name = "Size") +
            scale_y_continuous(breaks=seq(0,1,0.1)) +
            scale_x_continuous(breaks=seq(-10,50,10), limits=c(-10,50)) +
            plot_theme + xlab("Change") + ylab("Power") + ggtitle ("BiSeq")
            ggsave(paste0('%(base)s',"_change_vs_size_power_lineplot.png"),
            plot = p, width = 5, height = 5)

            agg_DMR_count = aggregate(cat_df$DMR,by=list(cat_df$size,
                cat_df$change), FUN=function(x) length(x))
            colnames(agg_DMR_count) = c("size","change","count")

            p = ggplot(agg_DMR_count, aes(y=as.numeric(change),
            x=as.numeric(as.character(size)), fill=as.numeric(count))) +
            geom_tile() +
            scale_fill_continuous(name = "Count") +
            scale_x_continuous(breaks=seq(1,11,2)) +
            scale_y_continuous(breaks=seq(-10,50,10), limits=c(-10,50)) +
            plot_theme + xlab("Size") + ylab("Change") + ggtitle ("BiSeq")
            ggsave(paste0('%(base)s',"_change_vs_size_count.png"),
            plot = p, width = 5, height = 5)
            sink()
            unlink("%(outfile)s")}''' % locals())
    BiSeq_power_analysis()  # run r code
    out = open(outfile, "w")
    out.write("plotted scatter at:%(base)s_scatter_means.png\n" % locals())
    out.write("plotted heatmaps: %(base)s_change_vs_size_power_heatmap.png\n"
              % locals())
    out.write("plotted heatmaps: %(base)s_change_vs_size_power_lineplot.png\n"
              % locals())
    out.write("plotted heatmaps: %(base)s_change_vs_size_count.png\n"
              % locals())
    out.close()


@cluster_runnable
def spikeInClustersPlotM3D(infile, outfile, groups):

    analysis_df = pd.read_csv(infile, sep="\t")
    cluster_characteristics = ["contig", "size", "s", "e", "initial", "change"]
    seqnames_split = analysis_df['seqnames'].apply(lambda x: x.split("_"))
    seqnames_split = [["".join((x[0], x[1])), x[2], x[3], x[4], x[5], x[6]]
                      if len(x) > 6 else x for x in seqnames_split]
    cluster_values_df = pd.DataFrame(data=seqnames_split,
                                     index=list(range(1, len(seqnames_split) + 1)))
    cluster_values_df.columns = cluster_characteristics
    concat_df = pd.concat([analysis_df, cluster_values_df], axis=1)
    groups = [x[0] for x in groups]

    within = copy.copy(cluster_characteristics)
    M3D_within_col = []
    for x in groups:
        regex = "\S%s\d_vs_\S%s\d" % (x, x)
        M3D_within_col.extend([col for col in concat_df.columns
                               if re.search(regex, col)])
    within.extend(M3D_within_col)

    between = copy.copy(cluster_characteristics)
    pairs = [x for x in itertools.combinations(groups, 2)]
    M3D_between_col = []
    for x in pairs:
        regex = "\S%s\d_vs_\S%s\d" % (x[0], x[1])
        M3D_between_col.extend([col for col in concat_df.columns
                                if re.search(regex, col)])
    between.extend(M3D_between_col)

    between_df = concat_df.ix[:, between]
    within_df = concat_df.ix[:, within]

    r_df_between = pandas2ri.py2ri(between_df)
    r_df_within = pandas2ri.py2ri(within_df)
    base = P.snip(infile, ".out")
    PowerPlot = r('''library(ggplot2)
                  library(reshape)
                  function(df_between, df_within){
                  cluster_characteristics = c("contig" ,"size", "s", "e",
                  "initial", "change")
                  m_between_df=melt(df_between,id.var=cluster_characteristics)
                  m_between_df['between'] = 1
                  m_between_df$change = gsub("m", "-", m_between_df$change)
                  m_within_df=melt(df_within,id.var=cluster_characteristics)
                  m_within_df['between'] = 0
                  m_cat_df = rbind(m_between_df,m_within_df)
                  quant = (quantile(m_within_df$value,probs=0.95))
                  limited = m_within_df[m_within_df$value>quant,]

                  permute_func = function(x, len_x, vector, n){
                  rand = sample(vector, n, replace=TRUE)
                  return(length(rand[rand>x])/n)}

                  m_within_df$p_value = sapply(m_within_df$value,
                  function(x) permute_func(x, length(m_within_df$value),
                  m_within_df$value, 10000))

                  m_between_df$p_value=sapply(m_between_df$value,
                  function(x) permute_func(x, length(m_within_df$value),
                  m_within_df$value, 10000))
                  m_between_df$p_value_adj = p.adjust(
                  m_between_df$p_value , method ="BH")

                  agg_between = aggregate(m_between_df$p_value_adj,
                  by=list(m_between_df$size, m_between_df$change),
                  FUN=function(x) length(x[x<0.05])/length(x))
                  colnames(agg_between) = c("size","change","power")

                  text_theme = element_text(size=20)
                  text_theme_s = element_text(size=15)
                  text_theme_l = element_text(size=25)
                  plot_theme = theme(axis.title.x = text_theme,
                  axis.title.y = text_theme, axis.text.x = text_theme_s,
                  axis.text.y = text_theme_s, legend.text = text_theme_s,
                  legend.title = text_theme_s, title = text_theme_l)

                  p = ggplot(agg_between,aes(y=as.numeric(change),
                  x=as.numeric(as.character(size)),fill=as.numeric(power))) +
                  geom_tile() +
                  scale_fill_continuous(limits=c(0,1), name = "Power",
                                        low="grey95", high="mediumpurple4") +
                  scale_x_continuous(breaks=seq(1,11,2)) +
                  scale_y_continuous(breaks=seq(-20,50,10), limits=c(-20,50)) +
                  theme_bw() +
                  plot_theme + xlab("Size") + ylab("Change") + ggtitle ("M3D")
                  ggsave(paste0('%(base)s',"_change_vs_size_power_heatmap.png"),
                  plot = p, width = 5, height = 5)

                  p = ggplot(agg_between,aes(x=as.numeric(change),
                  y=as.numeric(power),
                  col=factor(as.numeric(as.character(size))))) +
                  geom_line(size=2) +
                  scale_colour_discrete(name = "Size") +
                  scale_y_continuous(breaks=seq(0,1,0.1)) +
                  scale_x_continuous(breaks=seq(-20,50,10), limits=c(-20,50)) +
                  theme_bw() +
                  plot_theme + xlab("Change") + ylab("Power") + ggtitle ("M3D")
                  ggsave(paste0('%(base)s',"_change_vs_size_power_lineplot.png"),
                  plot = p, width = 5, height = 5)

                  agg_between_count= aggregate(m_between_df$p_value_adj,
                  by=list(m_between_df$size,m_between_df$change),length)
                  colnames(agg_between_count) = c("size","change","count")

                  p = ggplot(agg_between_count,aes(y=as.numeric(change),
                  x=as.numeric(as.character(size)), fill=as.numeric(count))) +
                  geom_tile() +
                  scale_fill_continuous(name = "Count") +
                  scale_x_continuous(breaks=seq(1,11,2)) +
                  theme_bw() +
                  scale_y_continuous(breaks=seq(-20,50,10), limits=c(-20,50)) +
                  plot_theme + xlab("Size") + ylab("Change") + ggtitle ("M3D")
                  ggsave(paste0('%(base)s',"_change_vs_size_count.png"),
                  plot = p, width = 5, height = 5)}''' % locals())

    PowerPlot(r_df_between, r_df_within)
    out = open(outfile, "w")
    out.write("plotted heatmaps: %(base)s_change_vs_size_power_heatmap.png\n"
              % locals())
    out.write("plotted heatmaps: %(base)s_change_vs_size_power_lineplot.png\n"
              % locals())
    out.write("plotted heatmaps: %(base)s_change_vs_size_count.png\n"
              % locals())
    out.close()


@cluster_runnable
def runBiSeq(infiles, outfile, sample_id):
    # use this when following callMethylationStatus(i.e full run)
    # cov_infiles = infiles

    cov_infiles = [x for x in infiles if sample_id in x]
    out = open(outfile, "w")
    base = P.snip(os.path.abspath(outfile), ".tsv")
    basenames = [os.path.basename(x) for x in cov_infiles]
    samples = [re.sub(r".fastq\S+", "", x) for x in basenames]
    groups = [x.split("-")[1] for x in samples]
    c = '","'.join(cov_infiles)
    s = '","'.join(samples)
    g = '","'.join(groups)
    out.write(
        "basenames: %(basenames)s\ngroups: %(groups)s\nsamples: %(samples)s\n"
        % locals())
    out.write("c: %(c)s\ns: %(s)s\ng: %(g)s " % locals())
    out.close()

    grdevices = importr('grDevices')
    grdevices.png(file=base + "_test_vario.png")

    rcode = r('''
    library("BiSeq")
    library("ggplot2")
    function(){
    raw = readBismark(files = c("%(c)s"),DataFrame(group = c("%(g)s"),
    row.names = c("%(s)s")))
    colData(raw)$group <- factor(colData(raw)$group)
    rrbs_rel <- rawToRel(raw)
    rrbs.clust.unlim <- clusterSites(
    object = raw,groups = colData(raw)$group,perc.samples = 1,
    min.sites = 5,max.dist = 50,mc.cores=8)
    ind.cov <- totalReads(rrbs.clust.unlim) > 0
    quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
    rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
    predictedMeth <- predictMeth(object = rrbs.clust.lim,h=25,mc.cores=8)
    temp_df=data.frame(y=methLevel(predictedMeth)[order(rowData(predictedMeth)),1])
    temp_df["x"]=methReads(rrbs.clust.lim)/((totalReads(rrbs.clust.lim))[,1])
    temp_df["cov"] = totalReads(rrbs.clust.lim)[,1]
    write.table(temp_df, file = paste0("%(base)s", "_biseq_df_1.tsv"),
    sep = "\t", quote = F)

    betaResults <- betaRegression(formula = ~group, link = "probit",
    object = predictedMeth, type = "BR",mc.cores=8)
    predictedMethNull <- predictedMeth
    colData(predictedMethNull)$group.null <- rep(c(1,2), 3)
    betaResultsNull <- betaRegression(formula = ~group.null,
    link = "probit",object = predictedMethNull, type = "BR",mc.cores=8)
    vario <- makeVariogram(betaResultsNull)
    length_vario = length(vario$variogram[[1]][,2])
    estimated_sill = median(vario$variogram[[1]][-(1:(length_vario-50)),2])
    vario.sm <- smoothVariogram(vario, sill = estimated_sill)

    plot(vario$variogram[[1]],ylim=c(0,5))
    lines(vario.sm$variogram[,c("h", "v.sm")],col = "red", lwd = 1.5)

    vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)
    vario.sm$pValsList <- vario.aux$pValsList

    locCor <- estLocCor(vario.sm)

    # testClusters throws an error if no clusters found so need to catch error
    clusters.rej <- try(testClusters(locCor, FDR.cluster = 0.1), silent=TRUE)

    if (class(clusters.rej) == "try-error")
      {
      clusters.trimmed = data.frame(chr=NA, pos=NA, p.val=NA, meth.group1=NA,
                                    meth.group2=NA, meth.diff=NA, estimate=NA,
                                    std.error=NA, pseudo.R.sqrt=NA,
                                    cluster.id=NA, z.score=NA, pos.new=NA,
                                    p.li=NA)
      }

    else
      {
      clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = 0.1)
      DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)
      }

    # just in case we want to do something else with the results
    save.image()

    write.table(clusters.trimmed, paste0("%(base)s", "_trimmed_clusters.tsv"),
    sep="\t", quote=F, row.names=F)

    concat_vario = do.call(rbind, vario.sm$pValsList)
    write.table(concat_vario, paste0("%(base)s", "_vario_table_test.tsv"),
    sep="\t",quote = F,row.names = F)

    row_p = rowData(predictedMeth)
    df_ranges<- data.frame(seqnames=seqnames(row_p), starts=start(row_p),
    ends=end(row_p))
    ranges_pred_meth_df=cbind(df_ranges,methLevel(predictedMeth))

    write.table(ranges_pred_meth_df, paste0("%(base)s","_predicted_meth.tsv"),
    sep="\t", quote = F,row.names = F)}''' % locals())

    rcode()
    grdevices.dev_off()
