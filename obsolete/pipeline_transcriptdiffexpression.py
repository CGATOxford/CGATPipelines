##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline transcriptdiffexpression
===========================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

To do: update documentation and pipeline.ini. e.g supports sailfish too!

Overview
========
RNA-Seq differential expression analysis can, broadly speaking, be
performed at two levels. Gene-level and transcript-level.

As transcripts are the true unit of expression, differential
expression at the transcript-level is more ideal. However,
quantification of transcript-level expression is complicated by reads
which align to multiple transcripts from the same gene, especially
with short read techonologies.  In addition transcript-level
quantification may be hindered by inadequate genome annotation.

Kallisto and Salmon are transcript quantification tools which attempt
to quantify transcripts directly from the sequence reads by
lightweight alignment algorithms - referred to as
"pseduoaligning". This avoids the time-consuming step of aligning
genes to the reference genome but depends heavily on the quality of
the reference transcript geneset.

Kallisto and Salmon can bootstrap the transcript abundance
estimates. In order to identify differentially expression transcripts,
Sleuth uses these bootstrap values to estimate the transcript-wise
techincal variance which is subtracted from the total variance, thus
leaving an estimate of the remaining biological variance. Sleuth then
allows the user to fit a transcript-wise general linear model to the
expression data to identify transcripts which are signficantly
differentially expressed.

These tools require a reference transcript geneset. The easiest way to
generate this is to use the 'auto-generate' method (see pipeline.ini)
which uses the output of pipeline_annotations.py with user-defined
filtering. Alternatively, the geneset may be user-supplied (must be
called 'geneset.fa'). If you're not using the 'auto-generate' option
and you want to perform the simulation with a pre-mRNA fraction
included, you must also include a 'geneset_pre_mrna.fa' geneset with
pre-mRNA sequences.

To generate a geneset multi-fasta from a gtf, use the following:

       zcat geneset.gtf |
        awk '$3=="exon"'|
        cgat gff2fasta
        --is-gtf --genome-file=genome.fa --fold-at=60 -v 0
        --log=geneset.fa.log > geneset.fa;
        samtools faidx geneset.fa
        '''

To generate a geneset multi-fasta of pre-mRNAs from a gtf, use the following:

        zcat geneset.gtf |
        awk '$3 == "transcript"'|
        cgat gff2fasta
        --is-gtf --genome-file=genome.fa --fold-at 60 -v 0
        --log=geneset_pre_mrna.fa.log > geneset_pre_mrna.fa;
        samtools faidx geneset_pre_mrna.fa


Prior to the sample quantification, reads are simulated from the gene
set. This is a naive RNA-Seq simulation which does not simulate the
well known but viable biases from library preparation
sequencing. Reads are sampled uniformly at random across the
transcript model and sequencing errors introduced at random uniformly
across the reads, with the fragment length sampled from a user-defined
normal distribution. The main purpose of the simulation is to flag
transcripts which cannot be accurately quantified with "near-perfect"
RNA-Seq reads, although it may also be used to compare the accuracy of
the tools selected, with the caveat that the simulation does not model
real RNA-Seq samples well. The user should check the performance with
the simulated data to establish whether the geneset used is
suitable. For instance, it has been noted that inclusion of poorly
support transcripts leads to poorer quantification of well-supported
transcripts.

Note: If the transcripts.fa is not being generated within the
pipeline, you must ensure the suppled geneset.fa is sorted in gene_id
order (gtf2gtf --method=sort --sort-order=gene) and you must supply a
file mapping transcript ids to gene ids called transcript2gene.tsv
with the following form:

transcript_1   g_1
transcript_2   g_1
transcript_3   g_2
transcript_4   g_3


Principal targets
-----------------

simulation
    perform the simulation only

quantification
    compute all quantifications

full
    compute all quantifications and perform differential transcript
    expression testing


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_transcriptdiffexpression.py config

Input files
-----------

Sequence read input. Can be fastq or sra, single or paired end.

Design_files ("*.design.tsv") are used to specify sample variates. The
minimal design file is shown below, where include specifies if the
sample should be included in the analysis, group specifies the sample
group and pair specifies whether the sample is paired. Note, multiple
design files may be included, for example so that multiple models can
be fitted to different subsets of the data

(tab-seperated values)

sample    include    group    pair
WT-1-1    1    WT    0
WT-1-2    1    WT    0
Mutant-1-1    1    Mutant    0
Mutant-1-2    1    Mutant    0

If further variates need to be given, e.g the General linear model is
defined as ~ group + replicate, these can be specified in further columns:

sample    include    group    pair    replicate
WT-1-1    1    WT    0    1
WT-1-2    1    WT    0    2
Mutant-1-1    1    Mutant    0    1
Mutant-1-2    1    Mutant    0    2

For each design file, the pipeline.ini must specify a model and contrasts


Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* kallisto >= 0.42.1
* salmon >= 0.5.0
* sleuth >= 0.27.1

Pipeline output
===============

The main outputs of the pipeline are results tables and plots from the
differential expression analysis. Outputs are generated for each
*.design.tsv file and each contrast specfied and placed in DEresults.dir

  `results_[design]_counts.tsv"`
    counts table for all samples within the design

  `results_[design]_tpm.tsv"`
    Transcripts Per Million (tpm) table for all samples within the design

  `results_[design]_[contrast]_sleuth_ma.png`
    MA plot using sleuth function

DEresults.dir contains further plots summarising the sleuth analysis

  `results_[design]_[contrast]_sleuth_vars.png`
    technical vs. observed variance plot from sleuth

  `results_[design]_[contrast]_sleuth_mean_var.png`
    mean-variance plot from sleuth

  `results_[design]_[contrast]_MA_plot.png`
    MA plot from sleuth results table
    (for direction comparison with other methods)

  `results_[design]_[contrast]_volcano_plot.png`
    volcano plot from sleuth results table
    (for direction comparison with other methods)

The summary_plots directory contains further plots summarising the
expression estimates across the samples

# Mention Simulation results too!


Glossary
========

.. glossary::


Code
====

"""

# To do:

# Once Kallisto is upgraded > 0.42.2, include alignment stats from parsing sam

# Once sleuth is capable of performing within gene comparisons of
# transcript expression, add this analysis here too

# Add power test using counts2counts.py?

# add option to remove flagged transcripts from gene set

from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import sqlite3
import glob
import pandas as pd
import numpy as np
import itertools

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Expression as Expression

import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineTracks as PipelineTracks

import CGATPipelines.PipelineTranscriptDiffExpression as TranscriptDiffExpression

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))


# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.sra")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

Sample = PipelineTracks.AutoSample
DESIGNS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob("*.design.tsv"), "(\S+).design.tsv")

###############################################################################
# load designs
###############################################################################


@transform(["%s.design.tsv" % x.asFile() for x in DESIGNS],
           suffix(".tsv"),
           ".load")
def loadDesigns(infile, outfile):
    '''load design files into database'''
    # note group column needs renaming

    tmpfile = P.getTempFilename("/ifs/scratch")

    statement = "sed 's/group/_group/g' %(infile)s > %(tmpfile)s"
    P.run()

    P.load(tmpfile, outfile)
    os.unlink(tmpfile)

###############################################################################
# Create geneset
###############################################################################

if PARAMS["geneset_auto_generate"]:

    # TS: move to module file?
    @mkdir("index.dir")
    @originate("index.dir/transcript_ids.tsv")
    def identifyTranscripts(outfile):
        '''output a list of gene identifiers where biotype matches filter'''

        dbh = connect()

        table = os.path.basename(
            PARAMS["annotations_interface_table_transcript_info"])

        where_cmd = "WHERE (%s)" % " OR ".join(
            ["gene_biotype = '%s'" % x for x in
             PARAMS["geneset_gene_biotypes"].split(",")])

        if PARAMS["geneset_transcript_biotypes"]:

            t_biotypes = PARAMS["geneset_transcript_biotypes"].split(",")
            where_cmd += " AND (%s)" % " OR ".join(
                ["transcript_biotype = '%s'" % x for x in t_biotypes])

        if PARAMS["geneset_transcript_support"]:

            if PARAMS["geneset_random_removal"]:
                # if using random removal we don't want to filter on
                # transcript support here
                pass

            else:
                # TS: TSL is not given for all transcripts. Filtering here
                # will retain transcripts without TSL annotation

                # TS: I'm using substr because the tsl values also describe
                # the previous tsl and we only want the current tsl
                support_cmd = " OR ".join(
                    ["substr(transcript_support,1,4) = 'tsl%s'" % x
                     for x in range(1, PARAMS["geneset_transcript_support"] + 1)])

                # ensembl transcript support not given (e.g "NA") for
                # pseudogenes, single exon transcripts, HLA, T-cell
                # receptor, Ig transcripts.  Do we want to keep these in?
                na_support_cmd = "substr(transcript_support,1,5) = 'tslNA' "

                null_support_cmd = "transcript_support IS NULL"

                where_cmd += " AND (%s OR %s OR %s )" % (support_cmd,
                                                         na_support_cmd,
                                                         null_support_cmd)

        if PARAMS["geneset_random_removal"]:
            # TS:this section is for testing null random removal of transcripts
            # perform random removal of transcripts
            # remove the equivalent number as would be removed by
            # transcript support level filtering
            # note: in line with above, transcripts with tsl = NA are retained

            select_cmd = """ SELECT DISTINCT gene_id, transcript_id,
            transcript_support FROM annotations.%(table)s %(where_cmd)s
            """ % locals()

            select = dbh.execute(select_cmd)

            previous_gene = ""
            transcript_ids = []
            tsls = []

            # TS: remove these counts when this section has been checked...
            n_entries = 0
            n_transcripts = 0
            n_genes = 0
            n_low_support = 0
            n_high_support = 0
            tsl_NA = 0
            tsl_NULL = 0

            with IOTools.openFile(outfile, "w") as outf:
                outf.write("transcript_id\n")

                for entry in select:

                    n_entries += 1

                    gene_id, transcript_id, tsl = entry

                    if gene_id == previous_gene:
                        n_transcripts += 1

                        transcript_ids.append(transcript_id)

                        # some transcripts do not have a tsl
                        try:

                            tsls.append(
                                int(tsl.strip().split()[0].replace("tsl", "")))
                        except:
                            if tsl is None:
                                tsl_NULL += 1
                            else:
                                if tsl.strip().split()[0] == "tslNA":
                                    tsl_NA += 1

                    else:
                        count_below_threshold = len(
                            [x for x in tsls if x >
                             PARAMS["geneset_transcript_support"]])
                        count_above_threshold = len(
                            [x for x in tsls if x <=
                             PARAMS["geneset_transcript_support"]])
                        n_low_support += count_below_threshold
                        n_high_support += count_above_threshold
                        if count_below_threshold > 0:
                            # randomly remove transcripts

                            transcript_ids = np.random.choice(
                                transcript_ids,
                                size=len(transcript_ids) -
                                count_below_threshold,
                                replace=False)

                        # for some gene_ids, all transcripts may be removed!
                        if len(transcript_ids) > 0:
                            outf.write(
                                "\n".join((x for x in transcript_ids)) + "\n")

                        previous_gene = gene_id
                        transcript_ids = [transcript_id]
                        try:
                            tsls = [
                                int(tsl.strip().split()[0].replace("tsl", ""))]
                        except:
                            tsls = []
                            if tsl is None:
                                tsl_NULL += 1
                            else:
                                if tsl.strip().split()[0] == "tslNA":
                                    tsl_NA += 1

                        n_transcripts += 1
                        n_genes += 1

                # random select again for the last gene_id
                count_below_threshold = len(
                    [x for x in tsls if x >
                     PARAMS["geneset_transcript_support"]])
                count_above_threshold = len(
                    [x for x in tsls if x <=
                     PARAMS["geneset_transcript_support"]])
                n_low_support += count_below_threshold
                n_high_support += count_above_threshold

                if count_below_threshold > 0:

                    transcript_ids = np.random.choice(
                        transcript_ids,
                        size=len(transcript_ids) - count_below_threshold,
                        replace=False)

                if len(transcript_ids) > 0:
                    outf.write("\n".join((x for x in transcript_ids)))

            print(("# entries %i, # transcripts %i, # genes %i,"
                   "# low support %i, # high support %i,"
                   " # NA tsl %i, # NULL tsl %i" % (
                       n_entries, n_transcripts, n_genes, n_low_support,
                       n_high_support, tsl_NA, tsl_NULL)))

        else:
            # select all distinct transcript_ids passing filters
            select_cmd = """ SELECT DISTINCT transcript_id
            FROM annotations.%(table)s %(where_cmd)s""" % locals()

            select = dbh.execute(select_cmd)

            with IOTools.openFile(outfile, "w") as outf:
                outf.write("transcript_id\n")
                outf.write("\n".join((x[0] for x in select)) + "\n")

    @transform(identifyTranscripts,
               regex("index.dir/transcript_ids.tsv"),
               "index.dir/transcripts.gtf.gz")
    def buildGeneSet(mapfile, outfile):
        ''' build a gene set with only transcripts which pass filter '''

        geneset = PARAMS['annotations_interface_geneset_all_gtf']

        statement = '''
        zcat %(geneset)s
        | cgat gtf2gtf
        --method=filter
        --filter-method=transcript
        --map-tsv-file=%(mapfile)s
        --log=%(outfile)s.log
        | cgat gtf2gtf
        --method=sort
        --sort-order=gene+transcript
        --log=%(outfile)s.log
        | gzip
        > %(outfile)s
        '''
        P.run()

    @transform(buildGeneSet,
               suffix(".gtf.gz"),
               ".fa")
    def buildReferenceTranscriptome(infile, outfile):
        ''' build reference transcriptome from geneset'''

        genome_file = os.path.abspath(
            os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

        statement = '''
        zcat %(infile)s |
        awk '$3=="exon"'|
        cgat gff2fasta
        --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
        --log=%(outfile)s.log > %(outfile)s;
        samtools faidx %(outfile)s
        '''
        P.run()

    @transform(buildGeneSet,
               suffix(".gtf.gz"),
               ".pre_mRNA.fa")
    def buildReferencePreTranscriptome(infile, outfile):
        ''' build a reference transcriptome for pre-mRNAs'''

        if PARAMS['simulation_pre_mrna_fraction']:
            genome_file = os.path.abspath(
                os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

            statement = '''
            zcat %(infile)s |
            awk '$3 == "transcript"'|
            cgat gff2fasta
            --is-gtf --genome-file=%(genome_file)s --fold-at 60 -v 0
            --log=%(outfile)s.log > %(outfile)s;
            samtools faidx %(outfile)s
            '''
            P.run()

        else:
            P.touch(outfile)


else:
    # if a reference gtf is provided, just soft link to this
    assert os.path.exists("geneset.fa") > 0, (
        "if not auto generating a geneset, you must"
        "provide a geneset in a geneset.fa file")

    @mkdir("index.dir")
    @files("geneset.fa", "index.dir/transcripts.fa")
    def buildReferenceTranscriptome(infile, outfile):
        ''' link to the geneset provided'''
        P.clone(os.path.abspath(infile), os.path.abspath(outfile))

    @mkdir("index.dir")
    @files("geneset.fa", "index.dir/transcripts.pre_mRNA.fa")
    def buildReferencePreTranscriptome(infile, outfile):
        ''' build a reference transcriptome for pre-mRNAs'''
        if PARAMS['simulation_pre_mrna_fraction']:
            assert os.path.exists("geneset_pre_mRNA.fa") > 0, (
                "if not auto generating a geneset and simulating with"
                " a pre-mRNA fraction, you must provide a 'pre-mrna'"
                " geneset in a 'geneset_pre_mRNA.fa' file")
            P.clone(os.path.abspath(infile), os.path.abspath(outfile))

        else:
            P.touch(outfile)

###############################################################################
# build indexes
###############################################################################


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".kallisto.index")
def buildKallistoIndex(infile, outfile):
    ''' build a kallisto index'''

    job_memory = "12G"

    statement = '''
    kallisto index -i %(outfile)s -k %(kallisto_kmer)s %(infile)s
    '''

    P.run()


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    ''' build a salmon index'''

    job_memory = "2G"

    statement = '''
    salmon index %(salmon_index_options)s -t %(infile)s -i %(outfile)s
    -k %(salmon_kmer)s
    '''

    P.run()


@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".sailfish.index")
def buildSailfishIndex(infile, outfile):
    ''' build a sailfish index'''

    # sailfish indexing is more memory intensive than Salmon/Kallisto
    job_memory = "6G"

    statement = '''
    sailfish index --transcripts=%(infile)s --out=%(outfile)s
    --kmerSize=%(sailfish_kmer)s
    %(sailfish_index_options)s
    '''

    P.run()


@follows(mkdir("index.dir"),
         buildReferencePreTranscriptome,
         buildKallistoIndex,
         buildSalmonIndex,
         buildSailfishIndex)
def index():
    pass


###############################################################################
# Simulation
###############################################################################
# if not simulating, final task ('simulation') is empty
if PARAMS['simulation_run']:

    @mkdir("simulation.dir")
    @transform(buildReferenceTranscriptome,
               suffix(".fa"),
               "_transcript_kmers.tsv",
               output_dir="simulation.dir")
    def countTranscriptKmers(infile, outfile):
        ''' count the number of unique and non-unique kmers per transcript '''

        job_memory = PARAMS["simulation_kmer_memory"]

        statement = '''
        cgat fasta2unique_kmers
        --input-fasta=%(infile)s
        --method=transcript
        --kmer-size=%(kallisto_kmer)s
        -L %(outfile)s.log
        > %(outfile)s '''

        P.run()

    @transform(countTranscriptKmers,
               suffix(".tsv"),
               ".load")
    def loadTranscriptKmers(infile, outfile):
        ''' load the kmer counts'''

        options = "--add-index=id"
        P.load(infile, outfile, options=options)

    @mkdir("simulation.dir")
    @transform(buildReferenceTranscriptome,
               suffix(".fa"),
               "_gene_kmers.tsv",
               output_dir="simulation.dir")
    def countGeneKmers(infile, outfile):
        ''' count the number of unique and non-unique kmers per gene '''

        job_memory = PARAMS["simulation_kmer_memory"]

        if PARAMS["geneset_auto_generate"]:
            genemap = P.getTempFilename(shared=True)

            dbh = connect()
            select = dbh.execute('''
            SELECT DISTINCT transcript_id, gene_id FROM transcript_info''')

            with IOTools.openFile(genemap, "w") as outf:
                for line in select:
                    outf.write("%s\t%s\n" % (line[0], line[1]))

        else:
            assert os.path.exists("transcript2gene.tsv"), (
                "if you want to run the simulation on a user-supplied "
                "geneset, you need to supply a file mapping "
                "transcripts to genes " "called transcript2gene.tsv")
            genemap = "transcript2gene.tsv"

        statement = '''
        cgat fasta2unique_kmers
        --input-fasta=%(infile)s
        --method=gene
        --genemap=%(genemap)s
        --kmer-size=%(kallisto_kmer)s
        -L %(outfile)s.log > %(outfile)s '''

        P.run()

        if PARAMS["geneset_auto_generate"]:
            os.unlink(genemap)

    @transform(countGeneKmers,
               suffix(".tsv"),
               ".load")
    def loadGeneKmers(infile, outfile):
        ''' load the kmer counts'''

        options = "--add-index=id"
        P.load(infile, outfile, options=options)

    @mkdir("simulation.dir")
    @follows(buildReferenceTranscriptome,
             buildReferencePreTranscriptome)
    @files([(["index.dir/transcripts.fa",
              "index.dir/transcripts.pre_mRNA.fa"],
             ("simulation.dir/simulated_reads_%i.fastq.1.gz" % x,
              "simulation.dir/simulated_read_counts_%i.tsv" % x))
            for x in range(0, PARAMS["simulation_iterations"])])
    def simulateRNASeqReads(infiles, outfiles):
        ''' simulate RNA-Seq reads from the transcripts fasta file
        and transcripts pre-mRNA fasta file'''

        # TS: to do: add option to learn parameters from real RNA-Seq data
        # TS: move to module file. the statement is complicated by
        # neccesity for random order for some simulations
        infile, premrna_fasta = infiles
        outfile, outfile_counts = outfiles

        single_end_random_cmd = ""
        paired_end_random_cmd = ""

        if PARAMS["simulation_paired"]:
            outfile2 = outfile.replace(".1.gz", ".2.gz")
            options = '''
            --output-paired-end
            --output-fastq2=%(outfile2)s ''' % locals()

            if PARAMS["simulation_random"]:

                # need to randomised order but keep pairs in same position
                tmp_fastq1 = P.getTempFilename()
                tmp_fastq2 = P.getTempFilename()

                # randomise fastqs, gzip and replace
                paired_end_random_cmd = '''
                ; checkpoint ;
                paste <(zcat %(outfile)s) <(zcat %(outfile2)s) |
                paste - - - - | sort -R |
                awk -F'\\t' '{OFS="\\n"; print $1,$3,$5,$7 > "%(tmp_fastq1)s";
                print $2,$4,$6,$8 > "%(tmp_fastq2)s"}'; checkpoint ;
                rm -rf %(outfile)s %(outfile2)s; checkpoint;
                gzip -c %(tmp_fastq1)s > %(outfile)s; checkpoint;
                gzip -c %(tmp_fastq2)s > %(outfile2)s
                ''' % locals()

                os.unlink(tmp_fastq1)
                os.unlink(tmp_fastq2)

        else:
            options = ""

            if PARAMS["simulation_random"]:
                single_end_random_cmd = '''
                paste - - - - | sort -R | sed 's/\\t/\\n/g'| '''

        if PARAMS["simulation_random"]:
            # random shuffling requires all the reads to be held in memory!
            # should really estimate whether 4G will be enough
            job_memory = "4G"
        else:
            job_memory = "1G"

        job_threads = 2

        statement = '''
        cat %(infile)s |
        cgat fasta2fastq
        --premrna-fraction=%(simulation_pre_mrna_fraction)s
        --infile-premrna-fasta=%(premrna_fasta)s
        --output-read-length=%(simulation_read_length)s
        --insert-length-mean=%(simulation_insert_mean)s
        --insert-length-sd=%(simulation_insert_sd)s
        --counts-method=copies
        --counts-min=%(simulation_copies_min)s
        --counts-max=%(simulation_copies_max)s
        --sequence-error-phred=%(simulation_phred)s
        --output-counts=%(outfile_counts)s
        --output-quality-format=33 -L %(outfile)s.log
        %(options)s | %(single_end_random_cmd)s
        gzip > %(outfile)s %(paired_end_random_cmd)s'''

        P.run()

    @mkdir("simulation.dir/quant.dir/kallisto")
    @transform(simulateRNASeqReads,
               regex("simulation.dir/simulated_reads_(\d+).fastq.1.gz"),
               add_inputs(buildKallistoIndex),
               r"simulation.dir/quant.dir/kallisto/simulated_reads_\1/abundance.h5")
    def quantifyWithKallistoSimulation(infiles, outfile):
        ''' quantify trancript abundance from simulated reads with kallisto'''

        # TS more elegant way to parse infiles and index?
        infiles, index = infiles
        infile, counts = infiles

        # multithreading not supported until > v0.42.1
        job_threads = PARAMS["kallisto_threads"]
        job_memory = "8G"

        kallisto_options = PARAMS["kallisto_options"]

        if PARAMS["simulation_bootstrap"]:
            kallisto_bootstrap = PARAMS["kallisto_bootstrap"]
        else:
            kallisto_bootstrap = 0

        m = PipelineMapping.Kallisto()
        statement = m.build((infile,), outfile)

        P.run()

    @transform(quantifyWithKallistoSimulation,
               suffix(".h5"),
               ".tsv")
    def extractKallistoCountSimulation(infile, outfile):
        ''' run kalliso h5dump to extract txt file'''

        outfile_dir = os.path.dirname(os.path.abspath(infile))

        statement = '''kallisto h5dump -o %(outfile_dir)s %(infile)s'''

        P.run()

    @mkdir("simulation.dir/quant.dir/salmon")
    @transform(simulateRNASeqReads,
               regex("simulation.dir/simulated_reads_(\d+).fastq.1.gz"),
               add_inputs(buildSalmonIndex),
               r"simulation.dir/quant.dir/salmon/simulated_reads_\1/quant.sf")
    def quantifyWithSalmonSimulation(infiles, outfile):
        # TS more elegant way to parse infiles and index?
        infiles, index = infiles
        infile, counts = infiles

        job_threads = PARAMS["salmon_threads"]
        job_memory = "8G"

        salmon_options = PARAMS["salmon_options"]

        if PARAMS["salmon_bias_correct"]:
            salmon_options += " --biascorrect"

        salmon_libtype = "ISF"

        if PARAMS["simulation_bootstrap"]:
            salmon_bootstrap = PARAMS["salmon_bootstrap"]
        else:
            salmon_bootstrap = 0

        m = PipelineMapping.Salmon(PARAMS["salmon_bias_correct"])
        statement = m.build((infile,), outfile)

        P.run()

    @mkdir("simulation.dir/quant.dir/sailfish")
    @transform(simulateRNASeqReads,
               regex("simulation.dir/simulated_reads_(\d+).fastq.1.gz"),
               add_inputs(buildSailfishIndex),
               r"simulation.dir/quant.dir/sailfish/simulated_reads_\1/quant.sf")
    def quantifyWithSailfishSimulation(infiles, outfile):
        # TS more elegant way to parse infiles and index?
        infiles, index = infiles
        infile, counts = infiles

        job_threads = PARAMS["sailfish_threads"]
        job_memory = "8G"

        sailfish_options = PARAMS["sailfish_options"]
        sailfish_libtype = "ISF"

        if PARAMS["simulation_bootstrap"]:
            sailfish_bootstrap = PARAMS["sailfish_bootstrap"]
        else:
            sailfish_bootstrap = 0

        m = PipelineMapping.Sailfish()
        statement = m.build((infile,), outfile)

        P.run()

    @transform(quantifyWithSalmonSimulation,
               regex("(\S+)/quant.sf"),
               r"\1/abundance.tsv")
    def extractSalmonCountSimulation(infile, outfile):
        ''' rename columns and remove comment to keep file format the same
        as kallisto'''

        # note: this expects column order to stay the same

        with IOTools.openFile(infile, "r") as inf:
            lines = inf.readlines()

            with IOTools.openFile(outfile, "w") as outf:
                outf.write("%s\n" % "\t".join(
                    ("target_id", "length", "tpm", "est_counts")))

                for line in lines:
                    if not line.startswith("# "):
                        outf.write(line)

    @transform(quantifyWithSailfishSimulation,
               regex("(\S+)/quant.sf"),
               r"\1/abundance.tsv")
    def extractSailfishCountSimulation(infile, outfile):
        ''' rename columns and remove comment to keep file format the same
        as kallisto'''

        # note: this expects column order to stay the same

        with IOTools.openFile(infile, "r") as inf:
            lines = inf.readlines()

            with IOTools.openFile(outfile, "w") as outf:
                outf.write("%s\n" % "\t".join(
                    ("target_id", "length", "tpm", "est_counts")))

                for line in lines:
                    if not line.startswith("# "):
                        outf.write(line)

    # define simulation targets
    SIMTARGETS = []

    mapToSimulationTargets = {'kallisto': (extractKallistoCountSimulation, ),
                              'salmon': (extractSalmonCountSimulation, ),
                              'sailfish': (extractSailfishCountSimulation, )}

    for x in P.asList(PARAMS["quantifiers"]):
        SIMTARGETS.extend(mapToSimulationTargets[x])

    @follows(*SIMTARGETS)
    def quantifySimulation():
        pass

    @transform(SIMTARGETS,
               regex(
                   "simulation.dir/quant.dir/(\S+)/simulated_reads_(\d+)/abundance.tsv"),
               r"simulation.dir/quant.dir/\1/simulated_reads_\2/results.tsv",
               r"simulation.dir/simulated_read_counts_\2.tsv")
    def mergeAbundanceCounts(infile, outfile, counts):
        ''' merge the abundance and simulation counts files for
        each simulation '''

        TranscriptDiffExpression.mergeAbundanceCounts(
            infile, outfile, counts, job_memory="2G", submit=True)

    @collate(mergeAbundanceCounts,
             regex("simulation.dir/quant.dir/(\S+)/simulated_reads_\d+/results.tsv"),
             r"simulation.dir/\1_simulation_results.tsv")
    def concatSimulationResults(infiles, outfile):
        ''' concatenate all simulation results '''

        df = pd.DataFrame()

        for inf in infiles:
            df_tmp = pd.read_table(inf, sep="\t")
            df = pd.concat([df, df_tmp], ignore_index=True)

        df.to_csv(outfile, sep="\t", index=False)

    @transform(concatSimulationResults,
               suffix("results.tsv"),
               add_inputs(countTranscriptKmers),
               "correlations.tsv")
    def calculateCorrelations(infiles, outfile):
        ''' calculate correlation across simulation iterations per transcript'''

        TranscriptDiffExpression.calculateCorrelations(
            infiles, outfile, PARAMS['simulation_bin_step'],
            job_memory="8G", submit=True)

    @transform(calculateCorrelations,
               suffix(".tsv"),
               ".load")
    def loadCorrelation(infile, outfile):
        ''' load the correlations data table'''

        options = "--add-index=id"
        P.load(infile, outfile, options=options)

    @transform(calculateCorrelations,
               regex("simulation.dir/(\S+)_simulation_correlations.tsv"),
               r"simulation.dir/\1_flagged_transcripts.tsv")
    def identifyLowConfidenceTranscript(infile, outfile):
        '''
        identify the transcripts which cannot be confidently quantified
        these fall into two categories:

        1. Transcripts whose with poor accuracy of estimated counts

           - transcripts with >2 absolute fold difference between the
             sum of ground truths and the sum of estimated counts are
             flagged

        2. Transcripts with poor correlation between estimated counts

           - spline fitted to relationship between correlation and kmer fraction.
             cut-off of 0.9 used to define minimum kmer fraction threshold.
             transcripts below threshold are flagged

        2. is not yet implemented. Currently the minimum kmer fraction is
        hardcoded as 0.03. Need to implement automatic threshold
        generation from data
        '''

        job_memory = "2G"

        TranscriptDiffExpression.identifyLowConfidenceTranscripts(
            infile, outfile, submit=True)

    @transform(identifyLowConfidenceTranscript,
               suffix(".tsv"),
               ".load")
    def loadLowConfidenceTranscripts(infile, outfile):
        ''' load the low confidence transcripts '''

        options = "--add-index=transcript_id"
        P.load(infile, outfile, options=options)

    @mkdir("simulation.dir")
    @follows(loadTranscriptKmers,
             loadGeneKmers,
             loadCorrelation,
             loadLowConfidenceTranscripts)
    def simulation():
        pass

else:
    @follows(mkdir("simulation.dir"))
    def simulation():
        pass

###############################################################################
# Remove flagged transcripts
###############################################################################

# Add task to optionally remove flagged transcripts

###############################################################################
# Estimate transcript abundance
###############################################################################

# enable multiple fastqs from the same sample to be analysed together
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    SEQUENCEFILES_REGEX = regex(
        r"%s/%s\.(fastq.1.gz|fastq.gz|sra)" % (
            DATADIR, PARAMS["merge_pattern_input"].strip()))
    # the last expression counts number of groups in pattern_input
    SEQUENCEFILES_KALLISTO_OUTPUT = r"quant.dir/kallisto/%s/abundance.h5" % (
        PARAMS["merge_pattern_output"].strip())
    SEQUENCEFILES_SALMON_OUTPUT = r"quant.dir/salmon/%s/quant.sf" % (
        PARAMS["merge_pattern_output"].strip())
    SEQUENCEFILES_SAILFISH_OUTPUT = r"quant.dir/sailfish/%s/quant.sf" % (
        PARAMS["merge_pattern_output"].strip())
else:
    SEQUENCEFILES_REGEX = regex(
        r".*/(\S+).(fastq.1.gz|fastq.gz|sra)")
    SEQUENCEFILES_KALLISTO_OUTPUT = r"quant.dir/kallisto/\1/abundance.h5"
    SEQUENCEFILES_SALMON_OUTPUT = r"quant.dir/salmon/\1/quant.sf"
    SEQUENCEFILES_SAILFISH_OUTPUT = r"quant.dir/sailfish/\1/quant.sf"


@mkdir("quant.dir/kallisto")
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildKallistoIndex),
         SEQUENCEFILES_KALLISTO_OUTPUT)
def quantifyWithKallisto(infiles, outfile):
    ''' quantify trancript abundance with kallisto'''

    # TS more elegant way to parse infiles and index?
    infile = [x[0] for x in infiles]
    index = infiles[0][1]

    # multithreading not supported until > v0.42.1
    job_threads = PARAMS["kallisto_threads"]
    job_memory = "6G"

    kallisto_options = PARAMS["kallisto_options"]
    bootstrap = PARAMS["kallisto_bootstrap"]

    m = PipelineMapping.Kallisto(pseudobam=PARAMS['kallisto_psuedobam'])
    statement = m.build(infile, outfile)

    P.run()


@mkdir("quant.dir/salmon")
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSalmonIndex),
         SEQUENCEFILES_SALMON_OUTPUT)
def quantifyWithSalmon(infiles, outfile):
    # TS more elegant way to parse infiles and index?
    infile = [x[0] for x in infiles]
    index = infiles[0][1]

    job_threads = PARAMS["salmon_threads"]
    job_memory = "6G"

    salmon_options = PARAMS["salmon_options"]
    salmon_libtype = PARAMS["salmon_libtype"]

    bootstrap = PARAMS["salmon_bootstrap"]

    m = PipelineMapping.Salmon()
    statement = m.build(infile, outfile)

    P.run()


@mkdir("quant.dir/sailfish")
@collate(SEQUENCEFILES,
         SEQUENCEFILES_REGEX,
         add_inputs(buildSailfishIndex),
         SEQUENCEFILES_SAILFISH_OUTPUT)
def quantifyWithSailfish(infiles, outfile):
    # TS more elegant way to parse infiles and index?
    infile = [x[0] for x in infiles]
    index = infiles[0][1]

    job_threads = PARAMS["sailfish_threads"]
    job_memory = "6G"

    sailfish_options = PARAMS["sailfish_options"]
    sailfish_libtype = PARAMS["sailfish_libtype"]

    bootstrap = PARAMS["sailfish_bootstrap"]

    m = PipelineMapping.Sailfish()
    statement = m.build(infile, outfile)

    P.run()


@transform(quantifyWithSalmon,
           regex("quant.dir/(\S+)/(\S+)/quant.sf"),
           r"quant.dir/\1/\2/abundance.h5")
def convertFromSalmon(infile, outfile):
    ''' convert the ssalmon output to Sleuth compatible h5 file'''
    TranscriptDiffExpression.convertFromFish(infile, outfile)


@transform(quantifyWithSailfish,
           regex("quant.dir/(\S+)/(\S+)/quant.sf"),
           r"quant.dir/\1/\2/abundance.h5")
def convertFromSailfish(infile, outfile):
    ''' convert the sailfish output to Sleuth compatible h5 file'''
    TranscriptDiffExpression.convertFromFish(infile, outfile)


@transform((quantifyWithKallisto),
           regex("quant.dir/(\S+)/(\S+)/abundance.h5"),
           r"quant.dir/\1/\2/quant.sf")
def convertToFish(infile, outfile):
    ''' convert the kallisto output to a flatfile, as per
    sailfish/salmon output'''

    outfile_dir = os.path.dirname(outfile)

    statement = '''
    kallisto h5dump -o %(outfile_dir)s %(infile)s;
    mv %(outfile_dir)s/abundance.tsv %(outfile)s;
    rm -rf %(outfile_dir)s/bs_abundance_*.tsv'''

    P.run()


# define quantifier targets
QUANTTARGETS = []

mapToQuantificationTargets = {'kallisto': (quantifyWithKallisto, convertToFish),
                              'salmon': (quantifyWithSalmon, convertFromSalmon,),
                              'sailfish': (quantifyWithSailfish, convertFromSailfish,)}

for x in P.asList(PARAMS["quantifiers"]):
    QUANTTARGETS.extend(mapToQuantificationTargets[x])


@follows(*QUANTTARGETS)
def quantify():
    pass

###############################################################################
# Differential isoform expression analysis
###############################################################################


# T.S move to Expression.py
def estimateSleuthMemory(bootstraps, samples, transcripts):
    ''' The memory usage of Sleuth is dependent upon the number of
    samples, transcripts and bootsraps.

    A rough estimate is:
    24 bytes * bootstraps * samples * transcripts
    (https://groups.google.com/forum/#!topic/kallisto-sleuth-users/mp064J-DRfI)

    TS: I've found this to be a serious underestimate so we use a
    more conservative estimate here with a default of 2G
    '''

    estimate = (48 * PARAMS["kallisto_bootstrap"] * samples * transcripts)

    job_memory = "%fG" % (max(2.0, (estimate / 1073741824.0)))

    return job_memory


# note: location of transcripts fasta is hardcoded here!
def generate_sleuth_parameters_on_the_fly():

    quantifiers = P.asList(PARAMS["quantifiers"])
    designs = [x.asFile() for x in DESIGNS]

    parameters = []

    for design, quantifier in itertools.product(designs, quantifiers):
        outfiles = [
            r"DEresults.dir/%s_%s_sleuth_results.tsv" % (design, quantifier),
            r"DEresults.dir/%s_%s_sleuth_counts.tsv" % (design, quantifier),
            r"DEresults.dir/%s_%s_sleuth_tpm.tsv" % (design, quantifier)]
        parameters.append(
            ("%s.design.tsv" % design, outfiles, quantifier,
             "index.dir/transcripts.fa"))

    for job_parameters in parameters:
        yield job_parameters


@follows(*QUANTTARGETS)
@mkdir("DEresults.dir")
@files(generate_sleuth_parameters_on_the_fly)
def runSleuth(design, outfiles, quantifier, transcripts):
    ''' run Sleuth to perform differential testing '''

    outfile, counts, tpm = outfiles

    Design = Expression.ExperimentalDesign(design)
    number_samples = sum(Design.table['include'])

    number_transcripts = 0
    with IOTools.openFile(transcripts, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                number_transcripts += 1

    job_memory = estimateSleuthMemory(
        PARAMS["%(quantifier)s_bootstrap" % locals()],
        number_samples, number_transcripts)

    design_id = P.snip(design, ".design.tsv")

    model = PARAMS["sleuth_model_%s" % design_id]
    contrasts = PARAMS["sleuth_contrasts_%s" % design_id]

    outfile_pattern = P.snip(outfile, ".tsv")

    statement = '''
    cgat counts2table
    --design-tsv-file=%(design)s
    --output-filename-pattern=%(outfile_pattern)s
    --log=%(outfile_pattern)s.log
    --method=sleuth
    --fdr=%(sleuth_fdr)s
    --model=%(model)s
    --contrasts=%(contrasts)s
    --sleuth-counts-dir=quant.dir/%(quantifier)s
    --outfile-sleuth-count=%(counts)s
    --outfile-sleuth-tpm=%(tpm)s
    >%(outfile)s
    '''

    P.run()

    if PARAMS['sleuth_genewise']:

        # gene-wise sleuth seems to be even more memory hungry!
        # Use 2 * transcript memory estimate
        job_memory = estimateSleuthMemory(
            PARAMS["%(quantifier)s_bootstrap" % locals()],
            2 * number_samples, number_transcripts)

        outfile_genes = outfile.replace(".tsv", "_genes.tsv")
        counts_genes = counts.replace(".tsv", "_genes.tsv")
        tpm_genes = tpm.replace(".tsv", "_genes.tsv")
        outfile_pattern_genes = P.snip(outfile_genes, ".tsv")

        statement = '''
        cgat counts2table
        --design-tsv-file=%(design)s
        --output-filename-pattern=%(outfile_pattern_genes)s
        --log=%(outfile_pattern_genes)s.log
        --method=sleuth
        --fdr=%(sleuth_fdr)s
        --model=%(model)s
        --contrasts=%(contrasts)s
        --sleuth-counts-dir=quant.dir/%(quantifier)s
        --outfile-sleuth-count=%(counts_genes)s
        --outfile-sleuth-tpm=%(tpm_genes)s
        --sleuth-genewise
        --gene-biomart=%(sleuth_gene_biomart)s
        >%(outfile_genes)s '''

        P.run()


@follows(*QUANTTARGETS)
@mkdir("DEresults.dir")
@collate(QUANTTARGETS,
         regex("quant.dir/(\S+)/(\S+)/quant.sf"),
         [r"DEresults.dir/all_counts_\1.tsv",
          r"DEresults.dir/all_tpm_\1.tsv"])
def mergeCounts(infiles, outfiles):
    ''' run Sleuth on all samples to generate a full counts/tpm table'''

    counts, tpm = outfiles

    df_counts = pd.DataFrame()
    df_tpm = pd.DataFrame()

    def addSampleColumn(df, infile, index, column, sample,
                        index_name=None):
        ''' add a single column called <sample> from the infile'''

        tmp_df = pd.read_table(
            infile, sep="\t", usecols=[index, column],
            index_col=index,  comment="#")

        tmp_df.columns = [sample]

        if index_name:
            tmp_df.index.name = index_name

        df = pd.merge(df, tmp_df, left_index=True, right_index=True,
                      how="outer")
        return df

    for infile in infiles:
        sample = os.path.basename(os.path.dirname(infile))

        if "kallisto" in infile:
            df_counts = addSampleColumn(
                df_counts, infile, "target_id", "est_counts", sample,
                index_name="transcript_id")
            df_tpm = addSampleColumn(
                df_tpm, infile, "target_id", "tpm", sample,
                index_name="transcript_id")

        elif "salmon" in infile or "sailfish" in infile:
            df_counts = addSampleColumn(
                df_counts, infile, "Name", "NumReads", sample,
                index_name="transcript_id")
            df_tpm = addSampleColumn(
                df_tpm, infile, "Name", "TPM", sample,
                index_name="transcript_id")

    df_counts.to_csv(counts, sep="\t", index=True)
    df_tpm.to_csv(tpm, sep="\t", index=True)


@transform(mergeCounts,
           regex("DEresults.dir/all_(\S+)_(\S+).tsv"),
           r"DEresults.dir/all_gene_expression_\1_\2.tsv")
def aggregateCounts(infiles, outfile):
    ''' aggregate counts across transcripts for the same gene_id '''

    for infile in infiles:
        outfile = infile.split("_")
        outfile[1] = "gene_expression_" + outfile[1]
        outfile = "_".join(outfile)

        statement = '''SELECT DISTINCT transcript_id, gene_id FROM
                       annotations.transcript_info'''
        transcript_info_df = pd.read_sql(statement, connect())
        transcript_info_df.set_index("transcript_id", inplace=True)

        df = pd.read_table(infile, sep="\t")
        df.set_index("transcript_id", inplace=True)
        df = df.join(transcript_info_df)
        df = pd.DataFrame(df.groupby("gene_id").apply(sum))
        df.drop("gene_id", 1, inplace=True)

        df.to_csv(outfile, sep="\t", index=True)


@product(aggregateCounts,
         formatter(
             "DEresults.dir/all_gene_expression_counts_(?P<QUANTIFIER>\S+).tsv"),
         ["%s.design.tsv" % x.asFile() for x in DESIGNS],
         formatter(".*/(?P<DESIGN>\S+).design.tsv$"),
         "DEresults.dir/{DESIGN[1][0]}_{QUANTIFIER[0][0]}_deseq2_DE_results.tsv")
def diffExpressionDESeq2(infiles, outfile):
    counts_inf, design_inf = infiles

    outfile_pattern = P.snip(outfile, ".tsv")

    design_id = P.snip(design_inf, ".design.tsv")
    model = PARAMS["deseq2_model_%s" % design_id]
    contrasts = PARAMS["deseq2_contrasts_%s" % design_id]

    outfile_pattern = P.snip(outfile, ".tsv")

    # counts need to be rounded to ints first
    tmp_counts = P.getTempFilename(shared=True)
    df = pd.read_table(counts_inf, sep="\t", index_col=0)
    df = df.applymap(lambda x: round(x, 0))
    df.to_csv(tmp_counts, sep="\t", index=True)

    ref_group = PARAMS["deseq2_ref_group_%s" % design_id]
    if ref_group:
        ref_group_cmd = "--reference-group=%s" % ref_group
    else:
        ref_group_cmd = ""

    if PARAMS["deseq2_ihw"]:
        ihw = "--use-ihw"
    else:
        ihw = ""

    statement = '''
    cgat counts2table
    --tags-tsv-file=%(tmp_counts)s
    --design-tsv-file=%(design_inf)s
    --output-filename-pattern=%(outfile_pattern)s
    --log=%(outfile_pattern)s.log
    --method=deseq2
    --fdr=%(deseq2_fdr)s
    --model=%(model)s
    --contrasts=%(contrasts)s
    %(ihw)s
    %(ref_group_cmd)s
    >%(outfile)s
    '''

    P.run()
    os.unlink(tmp_counts)


@transform(diffExpressionDESeq2,
           suffix(".tsv"),
           ".load")
def loaddiffExpressionDESeq2(infile, outfile):
    TranscriptDiffExpression.loadSleuthTableGenes(
        infile, outfile,
        PARAMS["annotations_interface_table_gene_info"],
        PARAMS["geneset_gene_biotypes"],
        PARAMS["database"],
        PARAMS["annotations_database"])


@transform(aggregateCounts,
           suffix(".tsv"),
           ".load")
def loadGeneWiseAggregates(infile, outfile):
    ''' load tables from aggregation of transcript-level estimate to
    gene-level estimates'''

    P.load(infile, outfile, options="add-index=gene_id")
    P.load(infile.replace("counts", "tpm"), outfile.replace("counts", "tpm"),
           options="add-index=gene_id")


@transform(mergeCounts,
           suffix(".tsv"),
           ".load")
def loadMergedCounts(infiles, outfile):
    ''' load tables from Sleuth collation of all estimates '''

    for infile in infiles:
        outfile = infile.replace(".tsv", ".load")

        TranscriptDiffExpression.loadSleuthTable(
            infile, outfile,
            PARAMS["annotations_interface_table_transcript_info"],
            PARAMS["geneset_gene_biotypes"],
            PARAMS["database"],
            PARAMS["annotations_database"])


@transform(runSleuth,
           regex("(\S+)_(\S+)_(counts|tpm).tsv"),
           r"\1_\2_\3.load")
def loadSleuthTables(infile, outfile):
    ''' load tables from Sleuth '''

    TranscriptDiffExpression.loadSleuthTable(
        infile, outfile,
        PARAMS["annotations_interface_table_transcript_info"],
        PARAMS["geneset_gene_biotypes"],
        PARAMS["database"],
        PARAMS["annotations_database"])


@transform(runSleuth,
           suffix("_results.tsv"),
           "_DEresults.load")
def loadSleuthResults(infiles, outfile):
    ''' load Sleuth results '''
    # List of three Sleuth outputs received as infiles due to generator
    # First file is results list
    infile = infiles[0]
    tmpfile = P.getTempFilename("/ifs/scratch")

    table = os.path.basename(
        PARAMS["annotations_interface_table_transcript_info"])

    if PARAMS["geneset_gene_biotypes"]:
        where_cmd = "WHERE " + " OR ".join(
            ["gene_biotype = '%s'" % x
             for x in PARAMS["geneset_gene_biotypes"].split(",")])
    else:
        where_cmd = ""

    select = """SELECT DISTINCT
    transcript_id, transcript_biotype, gene_id, gene_name
    FROM annotations.%(table)s
    %(where_cmd)s""" % locals()

    df1 = pd.read_table(infile, sep="\t")
    df1.set_index("test_id", drop=False, inplace=True)

    df2 = pd.read_sql(select, connect())
    df2.set_index("transcript_id", drop=False, inplace=True)

    df = df1.join(df2)
    df.to_csv(tmpfile, sep="\t", index=True)

    options = "--add-index=transcript_id"
    P.load(tmpfile, outfile, options=options)
    os.unlink(tmpfile)


@follows(runSleuth,
         loadSleuthTables,
         loadSleuthResults,
         loadMergedCounts,
         loadGeneWiseAggregates,
         loaddiffExpressionDESeq2)
def differentialExpression():
    pass


###############################################################################
# Expression summary plots
###############################################################################

@mkdir("summary_plots")
@collate(runSleuth,
         regex("DEresults.dir/(\S+)_(\S+)_sleuth_results.tsv"),
         add_inputs(r"\1.design.tsv"),
         r"summary_plots/\1_\2_plots.log")
def expressionSummaryPlots(infiles, logfile):
    ''' make summary plots for expression values for each design file'''

    infiles, design_inf = infiles[0]
    results_inf, counts_inf, tpm_inf = infiles

    job_memory = "4G"

    TranscriptDiffExpression.makeExpressionSummaryPlots(
        tpm_inf, design_inf, logfile, submit=True, job_memory=job_memory)


###############################################################################
# Generic pipeline tasks
###############################################################################


@follows(differentialExpression,
         expressionSummaryPlots,
         simulation,
         loadDesigns)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
