"""
=====================================================
Community analysis of metgenomic shotgun sequencing
=====================================================


Pipeline_metagenomecommunities.py takes as input a set of fastq
files from a shotgun sequencing experiment of environmental samples
and assesses community structure and function.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions
(:term:`experiment`) with one or more biological and/or technical
replicates (:term:`replicate`). A :term:`replicate` within each
:term:`experiment` is a :term:`track`.

Community profiling
--------------------

The pipeline uses various tools for assessing the abundance of taxa
in an environmental sample. To assess relative abundance of taxa in
a community, we use metaphlan as it is easy to use and because it performs
alignments against a reduced set of clade-specific marker genes
it is relativelyfast. Metaphlan is likely to perform better where the
samples are derived fromhuman (e.g gut) as the database is significantly
overrepresented for human derived taxa.

Where metaphlan attempts to estimate taxa relative abundances it does not
attemtp to assign every read to a taxa. An alternative method is to use
kraken. Kraken utilises a megablast-like approach in order to search for
exact sequence matches between reads and sequences in the kraken database
(taxonomy). Like metaphlan, kraken assumes that there is little sequence
divergence between the sequenced samples and the data in the database. In
cases where sequences are derived from environmental samples that have not been
sequenced before this will result in very few sequences being assigned to
a taxa.

A third approach is to use a sensitive alignment algorithm.
At the time of writing this, a new approach to perform
sensitive blastx-like alignments was developed by the Huson lab.
It is called Diamond and is 16000 times faster than blast -
a requirement for large datasets. The pipeline utilises diamond
in cases where it is expected that there is high divergence between
sequences derived from the sample and the NCBI nr database.
Following alignment with diamond, the pipeline will attempt to assign
each read to a taxa using the lcammaper tool (lowest common ancestor (LCA))
also developed by the Huson lab.


Functional profiling
---------------------

Whether DNA-seq or RNA-seq is used, functional profiling can be performed.
The functional profiling techniques used in the pipeline rely on a set of
non-redundant gene sequences (amino acids). Diamond is used to sensitively
align reads to the non-redundant database e.g. MetaRef or IGC.


Differential abundance estimations
-----------------------------------

To detect differences in abundance of taxa or genes, we utilise the
metagenomeSeq R package from bioconductor. This package utilises a
zero-inflated gaussian micture model to compensate for undersampling of
taxa/genes bewteen samples - which may cause overestimation of differences
due to differences in library size.

see http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineDocumenation`). To
start with, use the files supplied with the :ref:`Example` data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while
``replicate`` denotes the :term:`replicate` within an
:term:`experiment`. The ``suffix`` determines the file type.  The
following suffixes/file types are possible:

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be
   sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.


Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------+
|*Program*           |*Version*          |*Purpose*   |
+--------------------+-------------------+------------+
|diamond             |                   | Sensitive  |
|                    |                   | alignment  |
|                    |                   | algorithm  |
+--------------------+-------------------+------------+
|lcamapper           |                   |Community   |
|                    |                   |profiling + |
|                    |                   |fuctional   |
|                    |                   |profiling   |
+--------------------+-------------------+------------+
|metaphlan           |                   |taxanomic   |
|                    |                   | relative   |
|                    |                   | abundance  |
|                    |                   | estimator  |
+--------------------+-------------------+------------+
|kraken              |                   |megablast   |
|                    |                   |taxanomic   |
|                    |                   |assignment  |
|                    |                   |of reads    |
+--------------------+-------------------+------------+
|metagenomeSeq       |                   |differential|
|                    |                   |abundance   |
|                    |                   |tool        |
+--------------------+-------------------+------------+


Pipeline output
===============

TODO::

Additional outputs are stored in the database file :file:`csvdb`.

Glossary
========

.. glossary::

Code
====

"""

# load modules
from ruffus import *

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
# import CGATPipelines.PipelineMetagenomeAssembly as PipelineMetagenomeAssembly
import CGATPipelines.PipelineMetagenomeCommunities \
    as PipelineMetagenomeCommunities
import CGAT.FastaIterator as FastaIterator
import CGAT.Metaphlan as Metaphlan
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import pysam
import CGAT.Fastq as Fastq
import pandas
# import CGATPipelines.PipelineTracks as PipelineTracks

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P

P.getParameters(["pipeline.ini",
                 "%s/pipeline.ini" % os.path.splitext(__file__)[0], ])


PARAMS = P.PARAMS

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################


# collect fastq.gz tracks
# TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
#     glob.glob("*.fastq.gz"), "(\S+).fastq.gz") +\
#     PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
#         glob.glob("*.fastq.1.gz"), "(\S+).fastq.1.gz")

# ALL = PipelineTracks.Sample3()
# EXPERIMENTS =
# PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
# CONDITIONS =
# PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
# TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

###################################################################
# sequence files as input
###################################################################
SEQUENCEFILES = ("*.fastq.gz", "*.fastq.1.gz", "*.fasta.gz")
SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fastq.gz|fastq.1.gz|fasta.gz)")

###################################################################
# connecting to database
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''
    dbh = sqlite3.connect(PARAMS["database"])
    return dbh

###################################################################
###################################################################
###################################################################
# load number of reads
###################################################################
###################################################################
###################################################################


@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"\1.nreads")
def countReads(infile, outfile):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()


@merge(countReads, "reads_summary.load")
def loadReadCounts(infiles, outfile):
    '''load read counts into database.'''
    to_cluster = False
    outf = P.getTempFile()
    outf.write("track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile(infile).readlines()
        nreads = int(lines[0][:-1].split("\t")[1])
        outf.write("%s\t%i\n" % (track, nreads))
    outf.close()
    inname = outf.name

    tablename = P.toTable(outfile)
    statement = '''cgat csv2db -t %(tablename)s --log=%(outfile)s.log
                  < %(inname)s > %(outfile)s'''
    P.run()
    os.unlink(outf.name)

# read count target


@follows(loadReadCounts)
def count_reads():
    pass

###################################################################
###################################################################
###################################################################
# Preprocessing reads for community analysis
###################################################################
###################################################################
###################################################################


@follows(mkdir("fasta.dir"))
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"fasta.dir/\1.fa.gz")
def preprocessReads(infile, outfile):
    '''
    create merged fasta file for use with metaphlan
    '''
    # check for second read in the pair
    if infile.endswith(".fastq.gz"):
        E.info("converting fastq file to fasta file")
        statement = '''fastq-to-fasta.py %(infile)s 2> %(outfile)s.log
                       | gzip > %(outfile)s'''
        P.run()

    elif infile.endswith(".1.gz"):
        read2 = P.snip(infile, ".1.gz") + ".2.gz"
        assert os.path.exists(read2), "file does not exist %s" % read2

        log = infile.replace("fastq.", "")
        statement = '''cgat fastqs2fasta
                   -a %(infile)s
                   -b %(read2)s
                   --log=%(log)s.log
                   | gzip > %(outfile)s'''
        P.run()

###################################################################
###################################################################
###################################################################
# Estimate taxonomic relative abundances using metaphlan
###################################################################
###################################################################
###################################################################


@active_if("metaphlan" in PARAMS.get("classifiers"))
@follows(mkdir("metaphlan.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"metaphlan.dir/\1.bt2out.txt")
def mapReadsWithMetaphlan(infile, outfile):
    '''
    map reads first with metaphlan against the marker
    database - will reduce running time for subsequent
    steps to assess abundances etc
    NOTE: IF PAIRED END, FILES WILL RUN USING FIRST READ IN PAIR
    '''
    db = PARAMS.get("metaphlan_db")
    nproc = PARAMS.get("metaphlan_nproc")
    options = PARAMS.get("metaphlan_bowtie2_options")
    assert os.path.exists(
        PARAMS["metaphlan_db"] + ".1.bt2"), \
        """missing file %s: Are you sure you have
           the correct database for bowtie2?""" \
        % PARAMS["metaphlan_db"] + ".1.bt2"
    statement = '''zcat %(infile)s | metaphlan.py %(infile)s
                                     --input_type multifastq
                                     --mpa_pkl %(metaphlan_pkl)s
                                     --bowtie2db %(db)s
                                     --nproc %(nproc)s
                                     --bt2_ps %(options)s
                                     --no_map
                                     --bowtie2out %(outfile)s
                                     &> %(outfile)s.log'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(mapReadsWithMetaphlan,
           regex("(\S+)/(\S+).bt2out.txt"),
           r"metaphlan.dir/\2.readmap")
def buildMetaphlanReadmap(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via
    blastn searching
    '''
    statement = '''metaphlan.py -t reads_map
                   --input_type bowtie2out %(infile)s
                   | cgat metaphlan2table
                     -t read_map
                    --log=%(outfile)s.log
                    > %(outfile)s; checkpoint
                    ; sed -i 's/order/_order/g' %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildMetaphlanReadmap, suffix(".readmap"), ".readmap.load")
def loadMetaphlanReadmaps(infile, outfile):
    '''
    load the metaphlan read maps
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@merge(loadMetaphlanReadmaps, "metaphlan.dir/taxonomic.counts")
def countMetaphlanTaxonomicGroups(infiles, outfile):
    '''
    count the total number of species that
    were found by metaphlan
    '''
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\tcount\n")
    taxons = ["_order", "class", "family",
              "genus", "kingdom", "phylum", "species"]
    dbh = connect()
    cc = dbh.cursor()
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_readmap")
        for taxon in taxons:
            count = cc.execute(
                """SELECT COUNT(DISTINCT %s) FROM %s"""
                % (taxon, table)).fetchone()[0]
            outf.write("\t".join([track, taxon, str(count)]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@transform(mapReadsWithMetaphlan,
           regex("(\S+)/(\S+).bt2out.txt"),
           r"metaphlan.dir/\2.relab")
def buildMetaphlanRelativeAbundance(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via
    blastn searching
    '''
    statement = '''metaphlan.py -t rel_ab --input_type bowtie2out %(infile)s
                   | cgat metaphlan2table -t rel_ab
                    --log=%(outfile)s.log
                    > %(outfile)s; checkpoint
                    ; sed -i 's/order/_order/g' %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildMetaphlanRelativeAbundance, suffix(".relab"), ".relab.load")
def loadMetaphlanRelativeAbundances(infile, outfile):
    '''
    load the metaphlan relative abundances
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@merge(loadMetaphlanRelativeAbundances, "metaphlan.dir/taxonomic.abundances")
def buildMetaphlanTaxonomicAbundances(infiles, outfile):
    '''
    build a file that combines taxonomic abundances
    from each sample
    '''
    dbh = connect()
    cc = dbh.cursor()
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\ttaxon\tabundance\tidx\n")
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_relab")
        for data in cc.execute(
            """SELECT taxon_level,
                      taxon,
                      rel_abundance FROM %s""" % table).fetchall():
            idx = track.split("_")[1]
            outf.write(
                "\t".join([track, data[0], data[1], str(data[2]), idx]) + "\n")
    outf.close()

#########################################
# metaphlan target
#########################################


@active_if("metaphlan" in PARAMS.get("classifiers"))
@follows(loadMetaphlanRelativeAbundances,
         buildMetaphlanTaxonomicAbundances,
         countMetaphlanTaxonomicGroups,
         loadMetaphlanReadmaps)
def Metaphlan():
    pass

###################################################################
###################################################################
###################################################################
# Classify reads using Kraken and count taxonomic assignments
###################################################################
###################################################################
###################################################################


@follows(mkdir("kraken.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"kraken.dir/\1.classified.tsv.gz")
def classifyReadsWithKraken(infile, outfile):
    '''
    classify reads using kraken
    '''
    job_memory = "30G"
    kraken_db = PARAMS.get("kraken_db")
    temp = P.getTempFilename(".")
    statement = '''kraken --db %(kraken_db)s
                   --fastq-input
                   --gzip-compressed
                   %(infile)s > %(temp)s;
                   cat %(temp)s
                  | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(classifyReadsWithKraken,
           suffix(".classified.tsv.gz"),
           ".counts.tsv.gz")
def buildKrakenCounts(infile, outfile):
    '''
    build kraken counts table
    note that the output is produced using
    metaphlan2table but these are COUNTS and
    not relative abundance estimates. Therefore
    the file is passed through sed on the way out.
    '''
    kraken_db = PARAMS.get("kraken_db")
    temp = P.getTempFilename(".")
    statement = '''kraken-mpa-report
                   --db %(kraken_db)s
                   <(zcat %(infile)s)
                   > %(temp)s;
                   cat %(temp)s
                  | cgat metaphlan2table
                  -t rel_ab
                  --log=%(outfile)s.log
                  | sed 's/rel_abundance/count/g'
                  | gzip > %(outfile)s
                  ; rm -rf %(temp)s'''

    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildKrakenCounts, suffix(".tsv.gz"), ".kraken.load")
def loadKrakenCounts(infile, outfile):
    '''
    load kraken report
    '''
    P.load(infile, outfile, "--add-index=taxon")

###################################################################
###################################################################
###################################################################


@follows(mkdir("counts.dir"))
@split(loadKrakenCounts, "counts.dir/*.kraken.counts.tsv.gz")
def buildKrakenLevelCounts(infiles, outfiles):
    '''
    split counts by taxonomic levels
    '''
    for infile in infiles:
        tablename = P.toTable(os.path.basename(infile))
        track = P.snip(os.path.basename(infile), ".counts.kraken.load")
        levels = [
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"]

        dbh = connect()
        cc = dbh.cursor()

        for level in levels:
            outname = "counts.dir/" + \
                track + \
                ".%s.kraken.counts.tsv.gz" % level
            outf = IOTools.openFile(outname, "w")
            outf.write("taxa\tcount\n")
            for data in cc.execute("""SELECT taxon,
                                      count FROM %s
                                      WHERE taxon_level == '%s'
                                   """ % (tablename, level)).fetchall():
                outf.write("\t".join(map(str, data)) + "\n")
            outf.close()

###################################################################
###################################################################
###################################################################


@transform(buildKrakenLevelCounts, suffix(".tsv.gz"), ".load")
def loadKrakenLevelCounts(infile, outfile):
    '''
    load kraken counts
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@split(buildKrakenLevelCounts, "counts.dir/*kraken.aggregated.counts.tsv.gz")
def mergeKrakenCountsAcrossSamples(infiles, outfiles):
    '''
    merge counts into a single table across samples - input
    into metagenomeSeq
    '''
    levels = [
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"]
    for level in levels:
        prefixes = glob.glob(
            "counts.dir/*.%(level)s.kraken.counts.tsv.gz" % locals())
        prefixes = ",".join(
            [P.snip(
                os.path.basename(x),
                ".%(level)s.kraken.counts.tsv.gz" % locals()
            ) for x in prefixes])

        outname = os.path.join(
            "counts.dir", level + ".kraken.aggregated.counts.tsv.gz")

        statement = '''cgat combine_tables
                       --missing=0
                       --columns=1
                       --take=count
                       --glob=counts.dir/*.%(level)s.kraken.counts.tsv.gz
                       --prefixes=%(prefixes)s
                       --log=%(outname)s.log
                       | gzip > %(outname)s'''
        P.run()

###################################################################
###################################################################
###################################################################


@active_if(["kraken" in PARAMS.get("classifiers")])
@follows(mergeKrakenCountsAcrossSamples)
def Kraken():
    pass

###################################################################
###################################################################
###################################################################
# Assign reads to taxa using Lowest Common Ancestor (LCA).
# Initial alignment is done with Diamond
###################################################################
###################################################################
###################################################################


@follows(mkdir("diamond.dir"))
@transform(preprocessReads,
           regex("(\S+)/(\S+).fa.gz"),
           r"diamond.dir/\2.diamond.tsv.gz")
def runDiamondOnRawSequences(infile, outfile):
    '''
    diamond is an ultra fast equivalent to blastx. It takes
    fasta files as input
    At present it will run one sequence from paired files
    '''
    temp = P.getTempFilename(".")
    outtemp = P.getTempFilename(".")

    job_threads = PARAMS["diamond_threads"]
    job_memory = PARAMS["diamond_memory"]

    db = PARAMS["diamond_db"]
    diamond_options = PARAMS["diamond_options"]

    statement = '''zcat %(infile)s > %(temp)s.fastq;
                   checkpoint;
                   diamond blastx
                   --db %(db)s
                   --query %(temp)s.fastq
                   --daa %(temp)s.daa
                   --threads %(job_threads)s
                   --log
                   %(diamond_options)s
                   &> %(outfile)s.log;
                   diamond view -a %(temp)s.daa | gzip > %(outfile)s;
                   checkpoint;
                   rm -rf %(temp)s %(temp)s.fastq %(temp)s.daa
                '''
    P.run()


###################################################################
###################################################################
###################################################################


@transform(runDiamondOnRawSequences, suffix(".tsv.gz"), ".lca.gz")
def runLCA(infile, outfile):
    '''
    run the lowest common ancestor algorithm
    on the blast output to assign reads to
    taxa - from mtools. Runs with defaults at
    the moment.
    '''
    job_memory = "25G"

    # filtering options
    filter_list = P.asList(PARAMS.get("lca_filter"))
    if filter_list:
        filter_stmt = " | grep -v " + " | grep -v ".join(filter_list)
    else:
        filter_stmt = ""

    track = P.snip(outfile, ".lca.gz")
    gi2taxid = PARAMS.get("megan_gi2taxid")
    outf_tax = P.snip(outfile, ".gz")
    options = PARAMS.get("lca_options")
    statement = '''lcamapper.sh
                   -i %(infile)s
                   -f Detect
                    %(options)s
                   -gt %(gi2taxid)s
                   -o %(outf_tax)s > %(outfile)s.log
                   ; cat %(outf_tax)s
                  %(filter_stmt)s
                  | gzip > %(outfile)s
                  ; checkpoint
                  ; rm -rf %(outf_tax)s'''
    P.run()


###################################################################
###################################################################
###################################################################


@follows(mkdir("taxa_map.dir"))
@transform(runLCA, regex("(\S+)/(\S+).lca.gz"), r"taxa_map.dir/\2.map.gz")
def buildTaxaMap(infile, outfile):
    '''
    build a map from kingdom through species for
    each lca file - allows to map clades in downstream
    analysis
    '''
    statement = '''zcat %(infile)s
                   | cgat lca2table
                   --output-map
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@merge(buildTaxaMap, "taxa_map.dir/aggregated_taxa.map.gz")
def aggregateTaxaMaps(infiles, outfile):
    '''
    build a union of taxa mappings from kingdom through species
    '''
    found = []
    outf = IOTools.openFile(outfile, "w")
    outf.write("kingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for infile in infiles:
        inf = IOTools.openFile(infile)
        # skip header
        header = inf.readline()
        for line in inf.readlines():
            if line in found:
                continue
            else:
                found.append(line)
                outf.write(line)
    outf.close()

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".classified.gz")
def buildLCA(infile, outfile):
    '''
    tabulate LCA output into nice format. Per read assignment
    '''
    statement = '''zcat %(infile)s
                   | cgat lca2table
                     --summarise=individual
                     --log=%(outfile)s.log
                   | sed -e 's/order/_order/g'
                   | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".level.count")
def countLcaPerLevelTaxa(infile, outfile):
    '''
    count the number of taxa as found using LCA algorithm
    '''
    job_memory = "20G"
    statement = '''zcat %(infile)s |
                   cgat lca2table
                   --summarise=level-counts
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(countLcaPerLevelTaxa, suffix(".count"), ".count.load")
def loadCountLcaPerLevelTaxa(infile, outfile):
    '''
    load taxa level counts
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".counts.tsv.gz")
def buildLcaCounts(infile, outfile):
    '''
    count the number of taxa as found using LCA algorithm
    '''
    job_memory = "20G"
    statement = '''zcat %(infile)s |
                   cgat lca2table
                   --summarise=taxa-counts
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildLcaCounts,
           suffix(".diamond.counts.tsv.gz"),
           ".counts.diamond.load")
def loadLcaCounts(infile, outfile):
    '''
    load taxa level counts
    '''
    P.load(infile, outfile, "--add-index=taxa")

###################################################################
###################################################################
###################################################################


@follows(mkdir("counts.dir"))
@split(loadLcaCounts, "counts.dir/*.diamond.counts.tsv.gz")
def buildLcaLevelCounts(infiles, outfiles):
    '''
    split counts by taxonomic levels
    '''
    for infile in infiles:
        tablename = P.toTable(os.path.basename(infile))
        track = P.snip(os.path.basename(infile), ".counts.diamond.load")
        levels = [
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"]

        dbh = connect()
        cc = dbh.cursor()

        for level in levels:
            outname = "counts.dir/" + \
                track + ".%s.diamond.counts.tsv.gz" % level
            outf = IOTools.openFile(outname, "w")
            outf.write("taxa\tcount\n")
            for data in cc.execute(
                """SELECT taxa,
                          count
                          FROM %s
                          WHERE
                          level == '%s'""" % (tablename, level)).fetchall():
                outf.write("\t".join(map(str, data)) + "\n")
            outf.close()

###################################################################
###################################################################
###################################################################


@transform(buildLcaLevelCounts, suffix(".tsv.gz"), ".load")
def loadLcaLevelCounts(infile, outfile):
    '''
    load LCA per taxonomic level counts
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@split(buildLcaLevelCounts, "counts.dir/*.diamond.aggregated.counts.tsv.gz")
def mergeLcaCountsAcrossSamples(infiles, outfiles):
    '''
    merge counts into a single table across samples
    '''
    levels = [
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"]
    for level in levels:
        prefixes = glob.glob(
            "counts.dir/*.%(level)s.diamond.counts.tsv.gz" % locals())
        prefixes = ",".join(
            [P.snip(os.path.basename(x),
                    ".%(level)s.diamond.counts.tsv.gz"
                    % locals()) for x in prefixes])

        outname = os.path.join(
            "counts.dir", level + ".diamond.aggregated.counts.tsv.gz")

        statement = '''cgat combine_tables
                       --missing=0
                       --columns=1
                       --take=count
                       --glob=counts.dir/*.%(level)s.diamond.counts.tsv.gz
                       --prefixes=%(prefixes)s
                       --log=%(outname)s.log
                       | gzip > %(outname)s'''
        P.run()

###############################################
###############################################
###############################################

COUNT_DATA = []
classifiers = {"kraken": mergeKrakenCountsAcrossSamples,
               "lca": mergeLcaCountsAcrossSamples}
for classifier in P.asList(PARAMS.get("classifiers")):
    COUNT_DATA.append(classifiers[classifier])

###############################################
###############################################
###############################################


@jobs_limit(1, "R")
@transform(COUNT_DATA,
           suffix(".counts.tsv.gz"),
           ".counts.rarefied.tsv")
def rarefyTaxa(infile, outfile):
    '''
    rarefy to minimum counts
    '''
    PipelineMetagenomeCommunities.rarefy(infile,
                                         outfile,
                                         PARAMS.get("rarefy_sample"))

###############################################
###############################################
###############################################

rarefy = {0: (COUNT_DATA, ".counts.tsv.gz"),
          1: (rarefyTaxa, ".counts.rarefied.tsv")}
RAREFY_FUNC = rarefy[PARAMS.get("rarefy_rarefy_taxa")][0]
RAREFY_SUFFIX = rarefy[PARAMS.get("rarefy_rarefy_taxa")][1]


@jobs_limit(1, "R")
@transform(RAREFY_FUNC,
           suffix(RAREFY_SUFFIX),
           ".proportion.tsv")
def buildLcaProportionsAcrossSamples(infile, outfile):
    '''
    build the proportion of reads mapped to
    each taxoomic level per sample
    '''
    PipelineMetagenomeCommunities.buildLcaProportionsAcrossSamples(
        infile,
        outfile,
        dtype="taxa")

###############################################
###############################################
###############################################


@jobs_limit(1, "R")
@transform(buildLcaProportionsAcrossSamples,
           suffix(".tsv"),
           ".cumproportion.pdf")
def plotLcaProportionDistributions(infile, outfile):
    '''
    plot the cumulative proportions of taxa
    '''
    PipelineMetagenomeCommunities.plotProportionDistributions(infile,
                                                              outfile)

###############################################
###############################################
###############################################


@jobs_limit(1, "R")
@transform(buildLcaProportionsAcrossSamples,
           suffix(".tsv"),
           ".stackedbar.pdf")
def barchartLcaProportions(infile, outfile):
    '''
    barchart description of percent reads
    mapping to each taxon
    '''
    if PARAMS["heatmaps_order"]:
        order = PARAMS.get("heatmaps_order")
    else:
        order = open(infile).readline()[:-1].split("\t")[:-1]
        order.sort()
        order = ",".join(order)
    PipelineMetagenomeCommunities.barchartProportions(infile,
                                                      outfile,
                                                      order,
                                                      dtype="taxa")

###############################################
###############################################
###############################################


@follows(plotLcaProportionDistributions,
         barchartLcaProportions)
def proportions():
    pass

###############################################
###############################################
###############################################


@transform(mergeLcaCountsAcrossSamples, suffix(".tsv.gz"), ".load")
def loadAggregatedCounts(infile, outfile):
    '''
    load aggregated counts
    '''
    P.load(infile, outfile, "--add-index=taxa")

############################
# LCA target
############################


@follows(mergeLcaCountsAcrossSamples,
         aggregateTaxaMaps)
def Lca():
    pass

###################################################################
###################################################################
###################################################################
# Diversity analysis using the R package Vegan
###################################################################
###################################################################
###################################################################


COUNT_DATA = []
classifiers = {"kraken": mergeKrakenCountsAcrossSamples,
               "lca": mergeLcaCountsAcrossSamples}
for classifier in P.asList(PARAMS.get("classifiers")):
    COUNT_DATA.append(classifiers[classifier])


###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@follows(mkdir("diversity.dir"))
@transform(COUNT_DATA,
           regex("(\S+)/(\S+).counts.tsv.gz"),
           r"diversity.dir/\2.rarefaction.pdf")
def runRarefactionAnalysis(infile, outfile):
    '''
    run rarefaction analysis - sample to minimum sample count
    and calculate taxonomic richness
    '''
    f, to, step = PARAMS.get("rarefaction_from"), \
        PARAMS.get("rarefaction_to"), \
        PARAMS.get("rarefaction_step")
    rdir = PARAMS.get("rscriptsdir")
    PipelineMetagenomeCommunities.rarefactionCurve(infile,
                                                   outfile,
                                                   rdir,
                                                   f=f,
                                                   step=step)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(COUNT_DATA,
           regex("(\S+)/(\S+).counts.tsv.gz"),
           r"diversity.dir/\2.richness.sig")
def testRichness(infile, outfile):
    '''
    test significance of richness using kruskal wallis test
    '''
    rdir = PARAMS.get("rscriptsdir")
    sample = PARAMS.get("richness_sample")
    PipelineMetagenomeCommunities.testRichness(infile,
                                               outfile,
                                               rdir,
                                               sample)


###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@follows(mkdir("diversity.dir"))
@transform(COUNT_DATA,
           regex("(\S+)/(\S+).counts.tsv.gz"),
           r"diversity.dir/\2.diversity.tsv")
def buildDiversity(infile, outfile):
    '''
    build flat file with diversity calculation
    '''
    rdir = PARAMS.get("rscriptsdir")
    ind = PARAMS.get("diversity_index")
    PipelineMetagenomeCommunities.buildDiversity(infile,
                                                 outfile,
                                                 rdir,
                                                 ind=ind)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@follows(mkdir("diversity.dir"))
@transform(COUNT_DATA,
           regex("(\S+)/(\S+).counts.tsv.gz"),
           r"diversity.dir/\2.diversity.pdf")
def barplotDiversity(infile, outfile):
    '''
    barplot diversity between conditions
    '''
    rdir = PARAMS.get("rscriptsdir")
    ind = PARAMS.get("diversity_index")
    PipelineMetagenomeCommunities.barplotDiversity(infile,
                                                   outfile,
                                                   rdir,
                                                   ind=ind)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@follows(mkdir("diversity.dir"))
@transform(COUNT_DATA,
           regex("(\S+)/(\S+).counts.tsv.gz"),
           r"diversity.dir/\2.diversity.sig")
def testDiversity(infile, outfile):
    '''
    significance testing on community-wide diversity
    estimate
    '''
    rdir = PARAMS.get("rscriptsdir")
    ind = PARAMS.get("diversity_index")
    PipelineMetagenomeCommunities.testDiversity(infile,
                                                outfile,
                                                rdir,
                                                ind=ind)


@follows(testDiversity,
         testRichness,
         runRarefactionAnalysis,
         barplotDiversity,
         buildDiversity)
def diversity():
    pass

###################################################################
###################################################################
###################################################################
# Functional profiling - use diamond to align to non-redundant
# set of proteins
###################################################################
###################################################################
###################################################################


@follows(mkdir("genes.dir"))
@transform(preprocessReads,
           regex("(\S+)/(\S+).fa.gz"),
           r"genes.dir/\2.diamond.genes.tsv.gz")
def runDiamondOnGenes(infile, outfile):
    '''
    diamond is an ultra fast equivalent to blastx. It takes
    fastq files as input
    At present it will run one sequence from paired files
    '''
    temp = P.getTempFilename(".")

    job_threads = PARAMS["diamond_threads"]
    job_memory = PARAMS["diamond_memory"]

    db = PARAMS["genes_db"]
    diamond_options = PARAMS["genes_diamond_options"]

    statement = '''zcat %(infile)s > %(temp)s.fastq;
                   checkpoint;
                   diamond blastx
                   --db %(db)s
                   --query %(temp)s.fastq
                   --threads %(job_threads)s
                   --daa %(temp)s.daa
                   --log
                   %(diamond_options)s
                   &> %(outfile)s.log;
                   diamond view -a %(temp)s.daa | gzip > %(outfile)s;
                   checkpoint;
                   rm -rf %(temp)s %(temp)s.fastq %(temp)s.daa
                '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(runDiamondOnGenes, suffix(".tsv.gz"), ".counts.tsv.gz")
def buildDiamondGeneCounts(infile, outfile):
    '''
    build gene level counts
    '''
    job_memory = PARAMS.get("genes_memory")
    options = PARAMS.get("genes_count_options")
    statement = '''zcat %(infile)s |
                   cgat diamond2counts
                   %(options)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''

    P.run()


###################################################################
###################################################################
###################################################################

@transform(buildDiamondGeneCounts,
           suffix(".tsv.gz"),
           ".load")
def loadDiamondGeneCounts(infile, outfile):
    '''
    load gene counts
    '''
    P.load(infile, outfile, "--add-index=taxa")

###################################################################
###################################################################
###################################################################


@merge(buildDiamondGeneCounts, "genes.dir/gene_counts.tsv.gz")
def mergeDiamondGeneCounts(infiles, outfile):
    '''
    merge counts across sample datasets
    '''
    # USE THE SAME GLOB AS IN THE COMBINING TABLES SCRIPT
    # maintain correct order
    prefixes = [
        P.snip(os.path.basename(x), ".genes.counts.tsv.gz")
        for x in glob.glob("genes.dir/*.genes.counts.tsv.gz")]
    prefixes = ",".join(prefixes)

    statement = '''cgat combine_tables
                   --missing=0
                   --columns=1
                   --take=count
                   --glob=genes.dir/*.genes.counts.tsv.gz
                   --prefixes=%(prefixes)s
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()


###################################################################
###################################################################
###################################################################


@transform(mergeDiamondGeneCounts, suffix(".tsv.gz"), ".annotated.tsv.gz")
def annotatePathways(infile, outfile):
    '''
    annotate NOGs with there respecetive pathways (functional categories)
    '''
    PipelineMetagenomeCommunities.annotate(infile,
                                           outfile,
                                           PARAMS.get("pathways_geneset"))

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(annotatePathways,
           suffix(".annotated.tsv.gz"),
           ".proportion.pathways.tsv")
def buildPathwayProportionsAcrossSamples(infile, outfile):
    '''
    build the proportion of reads mapped to
    each taxoomic level per sample
    '''
    PipelineMetagenomeCommunities.buildLcaProportionsAcrossSamples(
        infile,
        outfile,
        dtype="pathway")


###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(mergeDiamondGeneCounts,
           suffix("_counts.tsv.gz"),
           ".proportion.tsv")
def buildGeneProportionsAcrossSamples(infile, outfile):
    '''
    build the proportion of reads mapped to
    each taxoomic level per sample
    '''
    PipelineMetagenomeCommunities.buildLcaProportionsAcrossSamples(
        infile,
        outfile,
        dtype="gene")

###############################################
###############################################
###############################################


@jobs_limit(1, "R")
@transform([buildGeneProportionsAcrossSamples,
            buildPathwayProportionsAcrossSamples],
           suffix(".tsv"),
           ".cumproportion.pdf")
def plotGeneProportionDistributions(infile, outfile):
    '''
    plot the cumulative proportions of taxa
    '''
    PipelineMetagenomeCommunities.plotProportionDistributions(infile,
                                                              outfile)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform([buildGeneProportionsAcrossSamples,
            buildPathwayProportionsAcrossSamples],
           suffix(".tsv"),
           ".stackedbar.pdf")
def barchartGeneProportions(infile, outfile):
    '''
    heatmap description of percent reads
    mapping to each taxon
    '''
    if PARAMS["heatmaps_order_genes"]:
        order = PARAMS.get("heatmaps_order_genes")
    else:
        order = open(infile).readline()[:-1].split("\t")[:-1]
        order.sort()
        order = ",".join(order)
    PipelineMetagenomeCommunities.barchartProportions(infile,
                                                      outfile,
                                                      order,
                                                      dtype="pathways")

###################################################################
###################################################################
###################################################################


#################
# genes target
#################


@follows(mergeDiamondGeneCounts,
         loadDiamondGeneCounts,
         barchartGeneProportions)
def Genes():
    pass

###################################################################
###################################################################
###################################################################
# Count alignments in for each of the methods
###################################################################
###################################################################
###################################################################


COUNT_TARGETS = []
classifiers = {"kraken": loadKrakenLevelCounts,
               "lca": loadLcaLevelCounts}
for classifier in P.asList(PARAMS.get("classifiers")):
    COUNT_TARGETS.append(classifiers[classifier])

###################################################################


@follows(mkdir("alignment_stats.dir"))
@transform(COUNT_TARGETS + [loadDiamondGeneCounts],
           regex("(\S+)/(\S+).counts.load"),
           add_inputs(loadReadCounts),
           r"alignment_stats.dir/\2.stats")
def countAlignments(infiles, outfile):
    '''
    count queries that have been aligned and have
    and assignment
    '''
    infile = infiles[0]
    summary_table = P.toTable(infiles[1])

    # assume that files are named without any other R[0-9]
    track = os.path.basename(re.match(
        "(.*-R[0-9]*).(.*.counts.load)", infile).groups()[0])
    table = P.toTable(infile)

    # connect to database
    dbh = connect()
    cc = dbh.cursor()

    alignments = cc.execute("""SELECT SUM(count)
                               FROM %(table)s""" % locals()).fetchone()[0]
    nreads = cc.execute("""SELECT total_reads
                           FROM %(summary_table)s
                           WHERE track == '%(track)s'
                        """ % locals()).fetchone()[0]

    outf = open(outfile, "w")
    outf.write("total_reads\taligned_reads\tpaligned_reads\n")
    outf.write("\t".join(
        map(str,
            [nreads, alignments, (float(alignments) / nreads) * 100])) + "\n")
    outf.close()


###################################################################
###################################################################
###################################################################


@merge(countAlignments, "alignment_stats.dir/alignment_stats.tsv")
def aggregateAlignmentStats(infiles, outfile):
    '''
    one table for alignment stats
    '''
    samples = [
        os.path.basename(x).split(".")[0] for x in infiles if "genes" not in x]
    levels = [
        os.path.basename(x).split(".")[1] for x in infiles if "genes" not in x]
    tools = [
        os.path.basename(x).split(".")[2] for x in infiles if "genes" not in x]
    outf = IOTools.openFile(outfile, "w")
    outf.write("sample\tlevel\ttool\tnreads\tnassigned\tpassigned\n")
    for s, l, t in zip(samples, levels, tools):
        inf = IOTools.openFile(
            "alignment_stats.dir/" + s + "." + l + "." + t + "." + "stats")
        inf.readline()
        result = inf.readline()
        outf.write("\t".join([s, l, t]) + "\t" + result)
    outf.close()

###################################################################
###################################################################
###################################################################


@transform(aggregateAlignmentStats, suffix(".tsv"), ".pdf")
def plotAlignmentStats(infile, outfile):
    '''
    barplot alignment stats
    '''
    PipelineMetagenomeCommunities.barplotAlignmentStats(infile,
                                                        outfile,
                                                        take="passigned")

###################################################################
###################################################################
###################################################################


@follows(count_reads)
@transform(aggregateAlignmentStats, suffix(".stats"), ".stats.load")
def loadAlignmentStats(infile, outfile):
    '''
    load alignment counts
    '''
    P.load(infile, outfile)


@follows(loadAlignmentStats, plotAlignmentStats)
def Alignment_stats():
    pass

###################################################################
###################################################################
###################################################################
# Differential abundance analysis of taxa and genes. We use
# metagenomeSeq and DESeq2 here to assess differential abundance
###################################################################
###################################################################
###################################################################


CLASSIFIER_TARGETS = []
classifiers = {"kraken": mergeKrakenCountsAcrossSamples,
               "lca": mergeLcaCountsAcrossSamples}
for classifier in P.asList(PARAMS.get("classifiers")):
    CLASSIFIER_TARGETS.append(classifiers[classifier])

###################################################################


@follows(mkdir("diff.dir"))
@transform(CLASSIFIER_TARGETS + [mergeDiamondGeneCounts],
           regex("(\S+)/(\S+).tsv.gz"),
           r"diff.dir/\2.diff.tsv")
def runMetagenomeSeq(infile, outfile):
    '''
    run metagenomeSeq - a tool for calculating significance
    based on gene counts
    '''
    rscriptsdir = PARAMS.get("rscriptsdir")
    rscript = PARAMS.get("metagenomeseq_rscript")
    prefix = P.snip(infile.replace("counts.dir", "diff.dir"), ".tsv.gz")

    if infile.find("gene") != -1:
        prefix = P.snip(infile.replace("genes.dir", "diff.dir"), ".tsv.gz")
        k = PARAMS.get("metagenomeseq_genes_k")
        a = PARAMS.get("metagenomeseq_genes_a")

        if PARAMS.get("metagenomeseq_genes_restrict"):
            restrict_file = PARAMS.get("metagenomeseq_genes_restrict_file")
            temp = P.getTempFile(".")
            genes = set([x[:-1] for x in open(restrict_file).readlines()])
            inf = IOTools.openFile(infile)
            header = inf.readline()
            temp.write(header)
            for line in IOTools.openFile(infile).readlines():
                data = line[:-1].split("\t")
                if data[0] in genes:
                    temp.write(line)
            temp.close()
            infile = temp.name
    else:
        k = PARAMS.get("metagenomeseq_taxa_k")
        a = PARAMS.get("metagenomeseq_taxa_a")

        # hack just to look at genus
        if PARAMS.get("metagenomeseq_taxa_restrict") and "genus" in infile:
            restrict_file = PARAMS.get("metagenomeseq_taxa_restrict_file")
            temp = P.getTempFile(".")
            taxa = set([x[:-1] for x in open(restrict_file).readlines()])
            inf = IOTools.openFile(infile)
            header = inf.readline()
            temp.write(header)
            for line in IOTools.openFile(infile).readlines():
                data = line[:-1].split("\t")
                if data[0] in taxa:
                    temp.write(line)
            temp.close()
            infile = temp.name

    statement = '''%(rscript)s %(rscriptsdir)s/run_metagenomeseq.R
                   -c %(infile)s
                   -p %(prefix)s
                   --k %(k)i
                   --a %(a)f > %(outfile)s.log'''

    P.run()

###################################################################
###################################################################
###################################################################


@follows(mkdir("deseq2.dir"))
@transform(CLASSIFIER_TARGETS + [mergeDiamondGeneCounts],
           regex("(\S+)/(\S+).tsv.gz"),
           r"deseq2.dir/\2.diff.tsv")
def runDESeq2(infile, outfile):
    '''
    run DESeq2 - a tool for calculating significance
    based on gene counts
    '''
    # build design as a temporary file
    design = P.getTempFile(".")
    design.write("track\tgroup\tinclude\tpair\n")
    samples = IOTools.openFile(infile).readline()[:-1].split("\t")
    samples = samples[1:]
    conditions = [x.split("-")[1] for x in samples]
    for i in range(len(samples)):
        design.write("%s\t%s\t1\t0\n" % (samples[i], conditions[i]))
    design.close()
    d = design.name

    # run DESeq2
    outpattern = P.snip(outfile, ".diff.tsv") + "_"
    fdr = PARAMS.get("deseq2_fdr")
    min_rowcounts = PARAMS.get("deseq2_filter_min_counts_per_row")
    min_samplecounts = PARAMS.get("deseq2_filter_min_counts_per_sample")
    percentile_rowsums = PARAMS.get("deseq2_filter_percentile_rowsums")
    statement = '''cgat runExpression
                   --method=deseq2
                   --outfile=%(outfile)s
                   --output-filename-pattern=%(outpattern)s
                   --fdr=%(fdr)s
                   --tags-tsv-file=%(infile)s
                   --design-tsv-file=%(d)s
                   --filter-min-counts-per-row=%(min_rowcounts)s
                   --filter-min-counts-per-sample=%(min_samplecounts)s
                   --filter-percentile-rowsums=%(percentile_rowsums)s
                   --log=%(outfile)s.log'''

    P.run()

    os.unlink(design.name)

###################################################################
###################################################################
###################################################################


@subdivide([runMetagenomeSeq, runDESeq2],
           regex("(\S+)/(\S+).diff.tsv*"),
           add_inputs(aggregateTaxaMaps),
           r"\1/*.\2.diff.tsv")
def splitResultsByKingdom(infiles, outfiles):
    '''
    split results by kingdom for downstream
    analysis
    TODO: Tidy this up and put into module with
    P.submit
    '''
    result, mapfile = infiles

    # metagenomeseq normalised values are in .norm.matrix
    # and deseq2 normalised values are in .rlog.tsv.gz
    if os.path.dirname(infiles[0]) == "deseq2.dir":
        matrix = P.snip(infiles[0], ".diff.tsv") + "_rlog.tsv.gz"
    else:
        matrix = P.snip(infiles[0], ".diff.tsv") + ".norm.matrix"

    hierarchy = PipelineMetagenomeCommunities.readHierarchy(mapfile)

    # need to do it for both the results
    # table and the normalised matrix file
    for inf in [result, matrix]:
        header = IOTools.openFile(inf).readline()
        for kingdom, taxa in hierarchy.items():
            if kingdom == "NA":
                kingdom = "other"
            else:
                kingdom = kingdom
            # sepcify new outfile name
            outf = os.path.join(
                os.path.dirname(result), kingdom + ".")
            outf = outf + os.path.basename(result)
            if inf.endswith(".matrix"):
                suffix = ".norm.matrix"
                # last column in the matrix file
                taxon_ind = -1
            elif inf.endswith("_rlog.tsv.gz"):
                suffix = "_rlog.tsv.gz"
                taxon_ind = 0
            else:
                suffix = None
                taxon_ind = 8
            if suffix:
                outf = outf.replace(".diff.tsv", suffix)
            outf = IOTools.openFile(outf, "w")
            outf.write(header)
            for line in IOTools.openFile(inf).readlines():
                data = line.strip("\n").split("\t")
                taxon = data[taxon_ind].replace('"', '')
                if taxon in taxa:
                    outf.write(line)
                else:
                    continue
            outf.close()

###################################################################
###################################################################
###################################################################


@transform([runMetagenomeSeq, runDESeq2, splitResultsByKingdom],
           suffix(".tsv"), ".load")
def loadDifferentialAbundance(infile, outfile):
    '''
    load differentially abundant features
    '''
    P.load(infile, outfile, "--allow-empty-file")

###################################################################
###################################################################
###################################################################


@follows(mkdir("annotations.dir"))
@transform(PARAMS.get("genes_annotation"),
           regex("(\S+)/(\S+).txt"),
           r"annotations.dir/\2.load")
def loadGeneAnnotations(infile, outfile):
    '''
    load annotations file
    '''
    P.load(infile, outfile, "--header=COG,description")

###################################################################
###################################################################
###################################################################


@follows(mkdir("pathways.dir"))
@transform(loadDifferentialAbundance,
           regex("genes.dir/(\S+).diff.load"),
           r"pathways.dir/foreground.tsv")
def buildForegroundGeneset(infile, outfile):
    '''
    build foreground data set for pathways analysis
    '''
    table = P.toTable(infile)
    dbh = connect()
    cc = dbh.cursor()
    result = {}
    groups = set()
    for group in cc.execute(
        """SELECT group1, group2
           FROM %(table)s""" % locals()).fetchall():
        groups.add(group)

    for group in groups:
        result[group[0] + "_vs_" + group[1]] = {}

    p_type = PARAMS.get("metagenomeseq_taxa_threshold_option")
    logfc = PARAMS.get("metagenomeseq_taxa_fc_threshold")
    p = PARAMS.get("metagenomeseq_taxa_p_threshold")
    if p_type == "p":
        p_type = "P_Value"
    elif p_type == "padj":
        p_type = "adj_P_Val"

    for group in groups:
        group1, group2 = group[0], group[1]
        for data in cc.execute(
            """SELECT taxa, %(p_type)s, logFC FROM %(table)s
               WHERE
               group1 == '%(group1)s'
               AND
               group2 == '%(group2)s'""" % locals()).fetchall():
            gene_id, pval, logFC = data
            if pval < p and abs(logFC) > logfc:
                result[group[0] + "_vs_" + group[1]][gene_id] = 1
            else:
                result[group[0] + "_vs_" + group[1]][gene_id] = 0

    df = pandas.DataFrame(result)
    df.to_csv(outfile, sep="\t", index_label="gene_id")


###################################################################
###################################################################
###################################################################


@transform(loadDifferentialAbundance,
           regex("genes.dir/(\S+).diff.load"),
           r"pathways.dir/background.tsv")
def buildBackgroundGeneset(infile, outfile):
    '''
    build background data set for pathways analysis
    '''
    table = P.toTable(infile)
    dbh = connect()
    cc = dbh.cursor()

    outf = open(outfile, "w")
    outf.write("gene_id\n")
    for data in cc.execute(
        """SELECT DISTINCT taxa
           FROM %(table)s""" % locals()).fetchall():
        outf.write(data[0] + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@split([buildForegroundGeneset,
        buildBackgroundGeneset,
        PARAMS.get("pathways_geneset")],
       "pathways.dir/*.overall")
def runPathwaysAnalysis(infiles, outfiles):
    '''
    run pathways analysis
    '''
    genes, background, gene2pathway = infiles
    statement = '''cgat runGO \
                   --background=%(background)s
                   --genes=%(genes)s \
                   --filename-input=%(gene2pathway)s \
                   -q BH \
                   --fdr \
                   --output-filename-pattern=\
                   pathways.dir/%%(set)s.%%(go)s.%%(section)s" \
                   > pathways.dir/pathways.log  \
                '''
    P.run()

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform([runMetagenomeSeq, runDESeq2, splitResultsByKingdom],
           suffix(".diff.tsv"),
           ".pca.tsv")
def runPCA(infile, outfile):
    '''
    run principle components analysis
    '''
    # metagenomeseq normalised values are in .norm.matrix
    # and deseq2 normalised values are in .rlog.tsv.gz
    if os.path.dirname(infile) == "deseq2.dir":
        inf = P.snip(infile, ".diff.tsv") + "_rlog.tsv.gz"
        rownames = 1
    else:
        inf = P.snip(infile, ".diff.tsv") + ".norm.matrix"
        rownames = len(open(inf).readline().strip("\n").split("\t"))
    if len(IOTools.openFile(inf).readlines()) <= 2:
        E.warn("Empty matrix %s: Check this is correct" % inf)
        P.touch(outfile)
        P.touch(outfile.replace(".tsv", ".ve.tsv"))
    else:
        PipelineMetagenomeCommunities.runPCA(inf, outfile, rownames=rownames)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(runPCA, suffix(".tsv"), ".pdf")
def plotPCA(infile, outfile):
    '''
    plot principle components
    '''
    # also outputted a separate file for variance explained
    scores, ve = infile, P.snip(infile, ".tsv") + ".ve.tsv"

    if os.path.getsize(scores) == 0:
        E.warn("Empty matrix %s: Check this is correct" % scores)
        P.touch(outfile)
    else:
        PipelineMetagenomeCommunities.plotPCA(scores,
                                              ve,
                                              outfile)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform([runMetagenomeSeq, runDESeq2], suffix(".diff.tsv"), ".mds.pdf")
def runMDS(infile, outfile):
    '''
    run MDS analysis
    '''
    # metagenomeseq normalised values are in .norm.matrix
    # and deseq2 normalised values are in _rlog.tsv.gz
    if os.path.dirname(infile) == "deseq2.dir":
        inf = P.snip(infile, ".diff.tsv") + "_rlog.tsv.gz"
    else:
        inf = P.snip(infile, ".diff.tsv") + ".norm.matrix"
    PipelineMetagenomeCommunities.plotMDS(inf, outfile)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform([runMetagenomeSeq, runDESeq2, splitResultsByKingdom],
           suffix(".diff.tsv"),
           ".mds.sig")
def runPermanova(infile, outfile):
    '''
    run permanova on euclidean distances
    '''
    # metagenomeseq normalised values are in .norm.matrix
    # and deseq2 normalised values are in _rlog.tsv.gz
    if os.path.dirname(infile) == "deseq2.dir":
        inf = P.snip(infile, ".diff.tsv") + "_rlog.tsv.gz"
        rownames = 1
    else:
        inf = P.snip(infile, ".diff.tsv") + ".norm.matrix"
        rownames = len(open(inf).readline().strip("\n").split("\t"))

    # only run if the file is not empty
    if len(IOTools.openFile(inf).readlines()) == 1:
        E.warn("Empty matrix %s: Check this is correct" % inf)
        P.touch(outfile)
    else:
        PipelineMetagenomeCommunities.testDistSignificance(inf,
                                                           outfile,
                                                           rownames=rownames,
                                                           method="adonis")

###################################################################
###################################################################
###################################################################


@follows(plotPCA,
         runPermanova)
def PCA():
    pass

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(mergeLcaCountsAcrossSamples,
           suffix(".tsv.gz"), ".barplot.png")
def barplotAbundances(infile, outfile):
    '''
    barplots abundances
    '''
    # the infile is a separate file output by
    # run_metagenomeseq = normalised counts
    PipelineMetagenomeCommunities.barplotAbundance(infile, outfile)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform(runMetagenomeSeq, suffix(".tsv"), ".ma.png")
def MAPlot(infile, outfile):
    '''
    ma plot the results
    '''
    if infile.find("gene") != -1:
        threshold_option = PARAMS.get("metagenomeseq_genes_threshold_option")
        p = PARAMS.get("metagenomeseq_genes_p_threshold")
        fc = PARAMS.get("metagenomeseq_genes_fc_threshold")

    else:
        threshold_option = PARAMS.get("metagenomeseq_taxa_threshold_option")
        p = PARAMS.get("metagenomeseq_taxa_p_threshold")
        fc = PARAMS.get("metagenomeseq_taxa_fc_threshold")

    # MAPlot for each group pair

    PipelineMetagenomeCommunities.MAPlot(infile,
                                         threshold_option,
                                         p,
                                         fc,
                                         outfile)

###################################################################
###################################################################
###################################################################


@jobs_limit(1, "R")
@transform([runMetagenomeSeq, splitResultsByKingdom],
           suffix(".tsv"),
           ".heatmap.pdf")
def plotDiffHeatmap(infile, outfile):
    '''
    plot differentially expressed genes on a heatmap
    '''
    norm_file = P.snip(infile, ".diff.tsv") + ".norm.matrix"

    if infile.find("gene") != -1:
        threshold_option = PARAMS.get("metagenomeseq_genes_threshold_option")
        p = PARAMS.get("metagenomeseq_genes_p_threshold")
        fc = PARAMS.get("metagenomeseq_genes_fc_threshold")

    else:
        threshold_option = PARAMS.get("metagenomeseq_taxa_threshold_option")
        p = PARAMS.get("metagenomeseq_taxa_p_threshold")
        fc = PARAMS.get("metagenomeseq_taxa_fc_threshold")

    PipelineMetagenomeCommunities.plotHeatmap(infile,
                                              norm_file,
                                              threshold_option,
                                              p,
                                              fc,
                                              outfile)

################################
# differential abundance target
################################


@follows(PCA, loadDifferentialAbundance)
def Differential_abundance():
    pass


###################################################################
###################################################################
###################################################################

# @transform(mergeDiamondGeneCounts, suffix(".tsv.gz"), ".diff.tsv")
# def runMetagenomeSeqOnGenes(infile, outfile):
#     '''
#     run metagenomeSeq - a tool for calculating significance
#     based on gene counts
#     '''
#     rscriptsdir = PARAMS.get("rscriptsdir")
#     rscript = PARAMS.get("metagenomeseq_rscript")
#     prefix = P.snip(infile, ".tsv.gz")

#     statement = '''%(rscript)s %(rscriptsdir)s/run_metagenomeseq.R
#                    -c %(infile)s
#                    -p %(prefix)s
#                    --k %(metagenomeseq_genes_k)i
#                    --a %(metagenomeseq_genes_a)i > %(outfile)s.log'''

#     P.run()

# ###################################################################
# ###################################################################
# ###################################################################

# @transform(runMetagenomeSeqOnGenes, suffix(".tsv"), ".annotated.tsv")
# def annotateDifferentialAbundance(infile, outfile):
#     '''
#     annotate differential abundance table with gene names etc
#     '''
#     annotation = PARAMS.get("genes_annotation")
#     PipelineMetagenomeCommunities.annotate(infile,
#                                            annotation,
#                                            outfile)

# ###################################################################
# ###################################################################
# ###################################################################

# @transform(annotateDifferentialAbundance, suffix(".tsv"), ".load")
# def loadDifferentialAbundance(infile, outfile):
#     '''
#     load the results of metagenomeSeq analysis
#     '''
#     P.load(infile, outfile, "--add-index=taxa")

# ###################################################################
# ###################################################################
# ###################################################################

# @transform(runMetagenomeSeqOnGenes, suffix(".diff.tsv"), ".mds.pdf")
# def runMDSOnGenes(infile, outfile):
#     '''
#     run MDS analysis on genes
#     '''
#     # the infile is a separate file output by
#     # run_metagenomeseq = normalised counts
#     inf = P.snip(infile, ".diff.tsv") + ".norm.matrix"
#     PipelineMetagenomeCommunities.plotMDS(inf, outfile)

# ###################################################################
# ###################################################################
# ###################################################################

# @transform(runMetagenomeSeqOnGenes, suffix(".tsv"), ".heatmap.png")
# def plotGenesDiffHeatmap(infile, outfile):
#     '''
#     plot differentially expressed genes on a heatmap
#     '''
#     norm_file = P.snip(infile, ".diff.tsv") + ".norm.matrix"

#     PipelineMetagenomeCommunities.plotHeatmap(infile,
#                                               norm_file,
#                                               PARAMS.get("metagenomeseq_p_threshold"),
#                                               PARAMS.get("metagenomeseq_fc_threshold"),
#                                               outfile)


# ###################################################################
# ###################################################################
# ###################################################################

# @transform(runMetagenomeSeqOnGenes, suffix(".tsv"), ".load")
# def loadMetagenomeSeqOnGenes(infile, outfile):
#     '''
#     laod differential abundance results
#     '''
#     P.load(infile, outfile, "--add-index=taxa")


# ###################################################################
# ###################################################################
# ###################################################################
# ## Counting KEGG associations
# ###################################################################
# ###################################################################
# ###################################################################
# @follows(mkdir("kegg.dir"))
# @transform(PARAMS.get("kegg_tre"),
#            regex("(\S+)/(\S+).tre"),
#            add_inputs(PARAMS.get("kegg_map")),
#            r"kegg.dir/\2.tsv")
# def buildKeggTable(infiles, outfile):
#     '''
#     build kegg table mapping KO identifiers to pathways. This is
#     based on the file that was downloaded with mtools (D.Huson)
#     '''
#     keggtre, keggmap = infiles
#     statement = '''cgat keggtre2table
#                    --kegg-tre=%(keggtre)s
#                    --map=%(keggmap)s
#                    --log=%(outfile)s.log
#                    > %(outfile)s'''
#     P.run()

# ###################################################################
# ###################################################################
# ###################################################################
# @transform(buildKeggTable, suffix(".tsv"), ".load")
# def loadKeggTable(infile, outfile):
#     '''
#     load KEGG table
#     '''
#     P.load(infile, outfile)

# ###################################################################
# ###################################################################
# ###################################################################
# @follows(mkdir("kegg.dir"))
# @transform(runLCA,
#            regex("(\S+)/(\S+).lca.gz"),
#            add_inputs(buildKeggTable),
#            r"kegg.dir/\2.kegg.counts")
# def countKeggAssociations(infiles, outfile):
#     '''
#     count the number of reads associted with Kegg pathways
#     '''
#     job_options = "-l mem_free=25G"
#     infile = infiles[0].replace(".lca", ".kegg")
#     kegg_table = infiles[1]
#     level = PARAMS.get("kegg_level")
#     statement = '''zcat %(infile)s |
#                    cgat lcakegg2counts
#                    --kegg-table=%(kegg_table)s
#                    --level=%(level)s
#                    --method=proportion
#                    --log=%(outfile)s.log
#                    > %(outfile)s'''
#     P.run()

# ###################################################################
# ###################################################################
# ###################################################################
# @transform(countKeggAssociations, suffix(".counts"), ".counts.load")
# def loadCountKeggAssociations(infile, outfile):
#     '''
#     load counts of KO associations
#     '''
#     P.load(infile, outfile, "--header=pathway,p_annotated_reads")

# ###################################################################
# ###################################################################
# ###################################################################
# @follows(mkdir("kegg.dir"))
# @transform(runLCA,
#            regex("(\S+)/(\S+).lca.gz"),
#            add_inputs(buildKeggTable),
#            r"kegg.dir/\2.kegg.ko.counts")
# def countKeggGenes(infiles, outfile):
#     '''
#     count the number of reads associated with KO identifiers
#     '''
#     job_options = "-l mem_free=25G"
#     infile = infiles[0].replace(".lca", ".kegg")
#     kegg_table = infiles[1]
#     level = PARAMS.get("kegg_level")
#     statement = '''zcat %(infile)s |
#                    cgat lcakegg2counts
#                    --kegg-table=%(kegg_table)s
#                    --level=D
#                    --method=proportion
#                    --log=%(outfile)s.log
#                    > %(outfile)s'''
#     P.run()

# ###################################################################
# ###################################################################
# ###################################################################
# @transform(countKeggGenes, suffix(".ko.counts"), ".ko.counts.load")
# def loadCountKeggGenes(infile, outfile):
#     '''
#     load counts of KO associations
#     '''
#     P.load(infile, outfile, "--header=KO,p_annotated_reads")

# #########################################
# # kegg target
# #########################################
# @follows(loadCountKeggAssociations)
# def kegg():
#     pass

#########################################
# full target
#########################################

@follows(Alignment_stats,
         Differential_abundance,
         diversity,
         plotPCA,
         proportions,
         Genes)
def full():
    pass

####################
# report building
####################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''
    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''
    E.info("updating documentation")
    P.run_report(clean=False)

#########################################
#########################################
#########################################


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
