"""
pipeline_peakcalling.py - Produce Peaklist from Bam Files
===================================================

:Author: Katy Brown, Charlie George & Adam Cribbs
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========

The aim of this pipeline is to create peaklists in :term:`bed` files from
aligned reads in :term:`bam` files that can then be taken on to downstream
analysis (e.g. motif identification, quantification of peaks etc.). Pipeline
also and generates QC  statistics that will inform you about the quality of
the peaksets generated.

Functionality
-------------

- Takes Paired-end or single end :term:`Bam` files you want to call peaks in
  (e.g. ChIP-Seq or ATAC-Seq samples and their appropriate 'input' controls).
- Runs peakcallers
- Runs ChIPQC R package for QC statistics
- Produces peak lists in bed files to takeforward for downstream analysis.


    Optional functions:
    -------------------
    - Filter :term:`Bam` files to remove:
            - Duplicates
            - Secondary alignments
            - Unpaired reads for paired-end files
            - Reads overlapping 'blacklisted' regions
            - Mapping quality (MAPQ) score
    - Pool input files for peakcalling to give better peakcalling results
      when inputs have poor coverage or lack sequening depth
    - Perform Irreproducible Discovery Rate (IDR) analysis (described
      further below) to get a consensus list of 'highly reproducible peaks'
      and assess replicate quaility.


NOTE: WARNINGS!!!!
------------------

1. IDR analysis may not be approprate for all type of peak file - It works
best with transcription factor ChIPs or methodologies producing 'narrow peaks'
or peaks with well defined boundaries.

'BroadPeak' IDR (e.g. for widespread histone marks such as H3K27ac)
might not work because peak boundary's are harder to define and thus may
not be so reproducible between replicates


2. Always check your output from this pipeline in a genome browser to check
peaks are being called suffiently!

3. This pipeline references :term:`ChIP bams` throughout in the code -this
references the immunoprecipitated (IP) sample from a ChIP experiment
(i.e. the file you want to find peaks in), :term:`Input bams` refer to the
bams of the input control samples that are used for background
normalisation in peak calling. Although we refer to :term:`ChIP bams`
this is only nomenclature and you could just as easily use
an ATAC-Seq :term:`bam` file or other :term:`bam` files in which you are
looking for peaks.

4) Whilst you can call peaks with as many peakcallers that are implemented in
the pipeline, only the results from one peakcaller can be taken forward for IDR
analysis. If you want to run IDR analysis on the output of multiple peakcallers
you will need first run IDR with one peakcaller then clone the pipeline, modify
pipeline.ini file and delete the appropriate files to rerun the IDR analysis on
the output from a different peakcaller. Be warned that IDR analysis generates
a large number of peakfiles and it's best to decide on your prefered peakcaller
before running the IDR analysis.


References
==========

This pipeline follows closely the ENCODE3 version 1 peakprocessing pipeline
described by Anshul Kundaje's group and the open source AQUAS TF ChIP-Seq
pipeline implemented by the Kundaje group:
    * (https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq)
    * (https://github.com/kundajelab/TF_chipseq_pipeline)

IDR analysis workflow is described here
    * (https://sites.google.com/site/anshulkundaje/projects/idr)

for troubleshooting/discussion of the IDR workflow see and extra documentation
see:
    * (https://groups.google.com/forum/#!forum/idr-discuss)

for ChIP and ATAC-Seq quality guidelines see:
    * (https://www.encodeproject.org/data-standards/)


IDR Analysis
============

IDR analysis is used to:
    * Give an indication of how reproducible the peaks that are produced by the
      peakcallers are within a single sample
    * Give an indication of how reproducible the peaks that are produced by the
      peakcallers are within biological replicates
    * produce a `conservative` peak list of highly reproducible peaks that
      can be taken forward to downstream analysis
    * produce an `oracle` peakset of the a large number of mostly reproducible
      peaks that can be taken forward to downstream analysis
    * sometimes the `conserative` and the `oracle` peakset will be the same
      list.
    * for further information on IDR analysis see the links above

Important notes:
IDR analysis requires peaks are called with a relaxed threshold to generate
a peaklist that contains (ideally) > 120,000 peaks that will contain
reproducible or 'true' peaks along with alot of irreproduible 'false' peaks.

Requirements
============

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+---------+------------+------------------------------------------------+
|*Program*|*Version*   |*Purpose*                                       |
+---------+------------+------------------------------------------------+
|samtools |>=0.1.16    |bam/sam file manipulation & stats               |
+---------+------------+------------------------------------------------+
|bedtools |            |working with intervals                          |
+---------+------------+------------------------------------------------+
|picard   |>=1.42      |duplication stats. The .jar files need to be in |
|         |            | your CLASSPATH environment variable.           |
+---------+------------+------------------------------------------------+
|macs2	  |>=2.1.1.    |peakcalling                                 	|
+---------+------------+------------------------------------------------+
|Conda	  |	           |		?????????????		                  	|
+---------+------------+------------------------------------------------+
|python   |>= 3.0      |run IDR analysis - currently set up in a        |
|         | 	       |conda enviroment that the pipeline calls	    |
+---------+------------+------------------------------------------------+
|IDR      |>= 2.0.2    |IDR analysis of peaks (bed files)               |
|         |            |from: (https://github.com/nboley/idr)           |
+---------+------------+------------------------------------------------+
|R        |            | used for QC stats                              |
+---------+------------+------------------------------------------------+
|ChIPQC   |            |                                                |
|R Package|            |                                                |
+---------+------------+------------------------------------------------+
|SICER
+---------+------------+------------------------------------------------+

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.


Pipeline Input
==============

Sample_bam = bam file you want to call peaks on (i.e. ChiP Bam or ATAC-Seq Bam)

Input_bam = control file used as background reference in peakcalling
(e.g. input file for ChIP-seq)

pipeline.ini = File containing paramaters and options for
running the pipeline

design.tsv = Design file based on design file for R package DiffBind
Has the following collumns:


Pipeline output
===============

The aim of this pipeline is to output a list of peaks that
can be used for further downstream analysis.

The pipeline generates several new directories containing
output files - these can roughly be grouped into XXX main
stages of the pipeline

1) filtered_bams.dir
   ---------------------
    Directory containing filtered bam files created by removing duplicates and
    filtering of origional bam files. These filtered bam files are then taken
    forward to IDR/peakcalling. If no filtering or deduplication is specified
    in the ini file then this directory will contain symbolic links to the
    origional bam files.

    Directory contains:
            * :term:`bams` files (and thier indexes) that have been filtered
            according to specifications in pipeline.ini
            * a number of log files detailing the number of reads that have been
            filtered out for each reason.
            * for paired-end samples a file with the frequency of fragment
              lengths (the distance between the paired reads 5' start positions)


2) IDR.dir
    -------
    Directory conatining the output files from IDR analysis
    IDR is currently only set up to use with macs2 because this
    is recomended by the authors of IDR. If you require IDR for broad
    peaks it is recomended to use the macs2 broad peaks setting.
    These include the lists of reproducible peaks and stats and
    QC tables summarising the output of the IDR analysis

    Directory contains:
            * IDR_inputs.dir
    This directory contains the files that are

    IDR_inputs.dir

    macs2.dir/

    peakcalling_bams.dir/

    peaks_for_IDR.dir/

    pooled_bams.dir/
    Peaksets:

    Conservative Peakset = Only obtained if IDR analysis run
    IDR analysis
    This analysis does a comparision on a pair of peak files to

    Tables
    Contained in the database are several tables used for QC and
    analysis




Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *
import sys
import os
import re
import csv
import sqlite3
import glob
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.PipelinePeakcallingScrum as PipelinePeakcalling
import CGAT.BamTools as Bamtools
import CGAT.Database as DB
import itertools
import pandas as pd
import numpy as np
import rpy2
from rpy2.robjects import r as R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


#########################################################################
# Load PARAMS Dictionary from Pipeline.innni file options ###############
#########################################################################
# load options from pipeline.ini file into PARAMS dictionary
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


PARAMS = P.PARAMS

# add parameters from annotations pipeline.ini
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    prefix="annotations_",
    update_interface=True))


# load IDR parameters into a dictionary to pass to the IDR step
# IDR requires multiple parameters from the PARAMS dictionary
idrPARAMS = dict()

# get IDR peakcaller and params
idrpc = PARAMS['peakcalling_idrpeakcaller']
idrPARAMS['idrsuffix'] = PARAMS["%s_idrsuffix" % idrpc]
idrPARAMS['idrcol'] = PARAMS["%s_idrcol" % idrpc]
idrPARAMS['idrcolname'] = PARAMS['%s_idrcolname' % idrpc]
idrPARAMS['useoracle'] = PARAMS['IDR_useoracle']


#######################################################################
# Check for design file & Match ChIP/ATAC-Seq Bams with Inputs ########
#######################################################################

# This section checks for the design table and generates:
# 1. A dictionary, inputD, linking each input file and each of the various
#    IDR subfiles to the appropriate input, as specified in the design table
# 2. A pandas dataframe, df, containing the information from the
#    design table.
# 3. INPUTBAMS: a list of control (input) bam files to use as background for
#    peakcalling.
# 4. CHIPBAMS: a list of experimental bam files on which to call peaks on.

# if design table is missing the input and chip bams  to empty list. This gets
# round the import tests

if os.path.exists("design.tsv"):
    df, inputD = PipelinePeakcalling.readDesignTable("design.tsv",
                                                     PARAMS['IDR_poolinputs'])
    INPUTBAMS = list(set(df['bamControl'].values))
    CHIPBAMS = list(set(df['bamReads'].values))

else:
    E.warn("design.tsv is not located within the folder")
    INPUTBAMS = []
    CHIPBAMS = []

# TODO we need to add code to pick up empty input and chipbams list and cause
# pipeline to throw an error


########################################################################
# Check if reads are paired end
########################################################################

if CHIPBAMS and Bamtools.isPaired(CHIPBAMS[0]) is True:
    PARAMS['paired_end'] = True
else:
    PARAMS['paired_end'] = False


#########################################################################
# Connect to database
#########################################################################

def connect():
    '''connect to database.
    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' %\
                (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


###########################################################################
# start of pipelined tasks
# 1) Preprocessing Steps - Filter bam files & generate bam stats
###########################################################################


@transform("design.tsv", suffix(".tsv"), ".load")
def loadDesignTable(infile, outfile):
    ''' load design.tsv to database '''
    P.load(infile, outfile)


@active_if(PARAMS['input'] != 0)
@follows(mkdir("filtered_bams.dir"))
@transform(INPUTBAMS, regex("(.*).bam"),
           [r"filtered_bams.dir/\1_filtered.bam",
            r"filtered_bams.dir/\1_counts.tsv"])
def filterInputBAMs(infile, outfiles):
    '''
    Applies various filters specified in the pipeline.ini to the bam file
    Currently implemented are filtering:
        unwanted contigs based on partial name matching
        unmapped reads
        unpaired reads
        duplicate reads
        secondary alignment reads
        reads below a mapping quality (MAPQ) score
        reads overlapping with blacklisted regions specified in bed file.
    '''
    filters = PARAMS['filters_bamfilters'].split(",")
    bedfiles = PARAMS['filters_bedfiles'].split(",")
    blthresh = PARAMS['filters_blacklistthresh']
    if blthresh != "":
        blthresh = float(blthresh)
    PipelinePeakcalling.filterBams(infile, outfiles, filters, bedfiles,
                                   blthresh,
                                   PARAMS['paired_end'],
                                   PARAMS['filters_strip'],
                                   PARAMS['filters_qual'],
                                   PARAMS['filters_contigs_to_remove'],
                                   PARAMS['filters_keepint'])


@follows(mkdir("filtered_bams.dir"))
@transform(CHIPBAMS, regex("(.*).bam"), [r"filtered_bams.dir/\1_filtered.bam",
                                         r"filtered_bams.dir/\1_counts.tsv"])
def filterChipBAMs(infile, outfiles):
    '''
    Applies various filters specified in the pipeline.ini to the bam file
    Currently implemented are filtering:
        unmapped reads
        unpaired reads
        duplicate reads
        secondary alignment reads
        reads below a mapping quality (MAPQ) score
        reads overlapping with blacklisted regions specified in bed file.
    '''
    filters = PARAMS['filters_bamfilters'].split(",")
    bedfiles = PARAMS['filters_bedfiles'].split(",")
    blthresh = PARAMS['filters_blacklistthresh']
    if blthresh != "":
        blthresh = float(blthresh)
    PipelinePeakcalling.filterBams(infile, outfiles, filters, bedfiles,
                                   float(blthresh),
                                   PARAMS['paired_end'],
                                   PARAMS['filters_strip'],
                                   PARAMS['filters_qual'],
                                   PARAMS['filters_contigs_to_remove'],
                                   PARAMS['filters_keepint'])


# ############################################################################
# ##### Filtering Stats and QC
# ############################################################################
@transform((filterChipBAMs, filterInputBAMs), suffix("_filtered.bam"),
           [r"\1_filtered.bam",
            r"\1_counts.tsv"])
def filteredBams(infiles, outfiles):
    ''' dummy task to collect filtered bams and counts.tsv tables
    for imput and chip file for downstream QC & Stats'''


@merge((filterChipBAMs, filterInputBAMs), "post_filtering_read_counts.tsv")
def mergeFilteringStats(infiles, outfile):
    '''
    Generates a table of read counts in each bam file after removal of:
    duplicates: duplicates reads
    secondary:  secondary alignment
    unpaired: unpaired reads
    unmapped: unmapped reads
    lowqual: low quality reads
    blacklist_xxx: reads in the blacklist file xxx
    contigs: removal of contigs that match patterns specified in ini file
    '''
    counts = [i[1] for i in infiles]
    bigtab = pd.DataFrame()
    for c in counts:
        tab = pd.read_csv(c, sep="\t")
        tab['Input_Bam'] = c.replace("_counts.tsv", ".bam").split("/")[-1]
        bigtab = bigtab.append(tab)
    bigtab = bigtab.rename(columns={'none': 'pre_filtering'})
    cs = []
    for c in bigtab.columns:
        if c.endswith(".bed"):
            c = "blacklist_%s" % c.split("/")[-1]
        cs.append(c)
    bigtab.columns = cs
    bigtab.to_csv(outfile, sep="\t", index=False)


@merge(mergeFilteringStats, "post_filtering_read_counts.load")
def loadFilteringStats(infile, outfile):
    '''load filtering stats to database'''
    P.load(infile, outfile)


@merge((filterChipBAMs, filterInputBAMs), "post_filtering_check.tsv")
def mergeFilteringChecks(infiles, outfile):
    '''take individual filering checks that detail the number of reads in the
    filtered bam file that are found for each flag that should have set in the
    filters and merge them to produce single table'''

    counts = [i[0].replace(".bam", ".filteringlog") for i in infiles]
    bigtab = pd.DataFrame()
    for c in counts:
        tab = pd.read_csv(c, sep="\t", index_col=0,  header=None)
        tab = tab.transpose()
        tab['Input_Filename'] = c.split("/")[-1].replace(".filteringlog",
                                                         "")
        bigtab = bigtab.append(tab)
    bigtab.to_csv(outfile, sep="\t", index=False)


@transform(mergeFilteringChecks, suffix(".tsv"), ".load")
def loadFilteringChecks(infile, outfile):
    '''load filtering stats to database '''
    P.load(infile, outfile)


@active_if(PARAMS['paired_end'])
@transform((filterChipBAMs, filterInputBAMs), suffix(".bam"),
           "_fraglengths.load")
def loadFragmentLengthDistributions(infiles, outfile):
    '''loads fragment length distributions into database - fragment length can
    only be computed if sample is paired-end if samples are not this function
    is not run'''
    infile = infiles[0].replace(".bam", ".fraglengths")
    if len(IOTools.openFile(infile).readlines()) != 0:
        P.load(infile, outfile)
    else:
        os.system("touch %s" % outfile)


@transform((filterChipBAMs, filterInputBAMs), suffix(".bam"),
           ".idxstats")
def getIdxstats(infiles, outfile):
    '''gets idxstats for bam file so number of reads per chromosome can
    be plotted later'''
    infile = infiles[0]
    statement = '''samtools idxstats %(infile)s > %(outfile)s''' % locals()
    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(getIdxstats, "idxstats_reads_per_chromosome.load")
def loadIdxstats(infiles, outfile):
    '''merge idxstats files into single dataframe and load
    to database

    Loads tables into the database
       * mapped_reads_per_chromosome

    Arguments
    ---------
    infiles : list
        list where each element is a string of the filename containing samtools
        idxstats output. Filename format is expected to be 'sample.idxstats'
    outfile : string
        Logfile. The table name will be derived from `outfile`.'''

    PipelineMappingQC.loadIdxstats(infiles, outfile)


@transform((filterChipBAMs, filterInputBAMs),
           suffix(".bam"),
           ".picard_stats")
def buildPicardStats(infiles, outfile):
    ''' build Picard alignment stats '''
    infile = infiles[0]
    reffile = os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa")

    PipelineMappingQC.buildPicardAlignmentStats(infile,
                                                outfile,
                                                reffile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardStats, "picard_stats.load")
def loadPicardStats(infiles, outfile):
    '''merge alignment stats into single tables.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)


@follows(loadFilteringStats,
         loadFilteringChecks,
         loadIdxstats,
         loadPicardStats)
def filtering():
    ''' dummy task to allow all the filtering of bams & collection of stats'''
    pass


# ### Make bigwigs of filtered bam files #####################################

@transform((filterChipBAMs, filterInputBAMs),
           suffix(".bam"),
           ".bw")
def buildBigWig(infile, outfile):
    '''build wiggle files from bam files.

    Generate :term:`bigWig` format file from :term:`bam` alignment file

    Parameters
    ----------
    infile : str
       Input filename in :term:`bam` format
    outfile : str
       Output filename in :term:`bigwig` format
    annotations_interface_contigs : str
       :term:`PARAMS`
       Input filename in :term:`bed` format

    '''
    inf = infile[0]
    # scale by Million reads mapped
    reads_mapped = Bamtools.getNumberOfAlignments(inf)
    scale = 1000000.0 / float(reads_mapped)
    tmpfile = P.getTempFilename()
    contig_sizes = PARAMS["annotations_interface_contigs"]
    job_memory = "3G"
    statement = '''bedtools genomecov
    -ibam %(inf)s
    -g %(contig_sizes)s
    -bg
    -scale %(scale)f
    > %(tmpfile)s;
    checkpoint;
    bedGraphToBigWig %(tmpfile)s %(contig_sizes)s %(outfile)s;
    checkpoint;
    rm -f %(tmpfile)s
    '''
    P.run()


###############################################################################
#
# 2) IDR  - preparation of files (pooled & pseudobams) for IDR
#
###############################################################################


# These steps are required for IDR and are only run if IDR is requested
if int(PARAMS['IDR_run']) == 1:
    @follows(mkdir("pooled_bams.dir"))
    @split(filterChipBAMs,
           r"pooled_bams.dir/*_pooled_filtered.bam")
    def makePooledBams(infiles, outfiles):
        '''
        IDR requires one bam file for each replicate and a pooled bam
        file of all replicates for a particular condition and tissue.
        This function generates the pooled bam files.
        '''
        cond_tissues = set(df['Condition'] + "_" + df['Tissue'])

        # Take each combination of tissues and conditions from the design
        # tables
        for ct in cond_tissues:
            p = ct.split("_")
            cond = p[0]
            tissue = p[1].split(".")[0]

            # identify and read all bam files for this combination of
            # tissue and condition
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]
            innames = subdf['bamReads'].values
            innames = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in innames])

            out = "pooled_bams.dir/%s_pooled_filtered.bam" % ct

            # Generate a merged, sorted, indexed bam file combining
            # all bam files for this tissue and condition
            PipelinePeakcalling.mergeSortIndex(innames, out)

    @active_if(PARAMS['IDR_poolinputs'] != "all" and PARAMS['input'] != 0)
    @follows(mkdir('IDR_inputs.dir'))
    @split(filterInputBAMs, "IDR_inputs.dir/*_pooled_filtered.bam")
    def makePooledInputs(infiles, outfiles):
        '''
        As pooled BAM files are used in the IDR, pooled input files also
        need to be generated - combined bam files of all the input bam
        files for this tissue.
        If you have chosen the "all" option for IDR_poolinputs in the
        pipeline.ini, this step is skipped, as all inputs are pooled for
        all IDR analyses.
        '''
        cond_tissues = set(df['Condition'] + "_" + df['Tissue'])

        # Take each combination of tissues and conditions from the design
        # tables
        for ct in cond_tissues:
            p = ct.split("_")
            cond = p[0]
            tissue = p[1].split(".")[0]
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]

            # find the inputs linked to any bam files for this combination of
            # tissues and conditions
            inputs = subdf['bamControl'].values
            inputs = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in inputs])
            out = "IDR_inputs.dir/%s_pooled_filtered.bam" % ct

            # generate a sorted, index, merged bam file for all of these
            # inputs
            PipelinePeakcalling.mergeSortIndex(inputs, out)

else:
    @transform(filterChipBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'filtered_bams.dir/\1.bam')
    def makePooledBams(infile, outfile):
        '''
        Dummy task if IDR not requested.
        '''
        pass

    @active_if(PARAMS['input'] != 0)
    @transform(filterInputBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'filtered_bams.dir/\1.bam')
    def makePooledInputs(infile, outfile):
        pass


if int(PARAMS['IDR_run']) == 1:
    @follows(mkdir("peakcalling_bams.dir"))
    @subdivide((filterChipBAMs, makePooledBams),
               regex("(.*)_bams.dir/(.*).bam"),
               [r"peakcalling_bams.dir/\2_pseudo_1.bam",
                r"peakcalling_bams.dir/\2_pseudo_2.bam",
                r"peakcalling_bams.dir/\2.bam"])
    def makePseudoBams(infiles, outfiles):
        '''
        Generates pseudo bam files each containing approximately 50% of reads
        from the original bam file for IDR self consistency analysis.
        Also generates a link to the original BAM file in the
        peakcalling_bams.dir directory.

        '''
        # makePooledBams generates a single output whereas filterChipBAMS
        # generates a bam file and a table - a list of outputs
        if isinstance(infiles, list):
            infile = infiles[0]
        else:
            infile = infiles

        pseudos = outfiles[0:2]
        orig = outfiles[2]

        PipelinePeakcalling.makeBamLink(infile, orig)

        PipelinePeakcalling.makePseudoBams(infile, pseudos,
                                           PARAMS['paired_end'],
                                           PARAMS['IDR_randomseed'],
                                           PARAMS['filters_bamfilters'].split(
                                               ","),
                                           submit=True)
else:
    @follows(mkdir('peakcalling_bams.dir'))
    @transform(filterChipBAMs, regex("filtered_bams.dir/(.*)_filtered.bam"),
               r'peakcalling_bams.dir/\1.bam')
    def makePseudoBams(infile, outfile):
        '''
        Link to original BAMs without generating pseudo bams
        if IDR not requested.
        '''
        PipelinePeakcalling.makeBamLink(infile[0], outfile)


# These three functions gather and parse the input (control) bam files into the
# IDR_inputs.dir directory prior to IDR analysis.
# The method used to do this depends on the IDR_poolinputs parameter

if PARAMS['IDR_poolinputs'] == "none":
    @active_if(PARAMS['input'] != 0)
    @follows(mkdir('IDR_inputs.dir'))
    @transform(filterInputBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'IDR_inputs.dir/\1.bam')
    def makeIDRInputBams(infile, outfile):
        '''
        When pooled inputs are not requested, the appropriate inputs are
        generated above in the filterInputBAMS step - this function links to
        these in the IDR_inputs.dir directory.
        '''
        infile = infile[0]
        PipelinePeakcalling.makeBamLink(infile, outfile)


elif PARAMS['IDR_poolinputs'] == "all":
    @active_if(PARAMS['input'] != 0)
    @follows(mkdir('IDR_inputs.dir'))
    @merge(filterInputBAMs, "IDR_inputs.dir/pooled_all.bam")
    def makeIDRInputBams(infiles, outfile):
        '''
        When all inputs are to be pooled and used as a control against all
        samples, a single merged bam is generated from the output of
        the filterInputBAMs step above in the IDR_inputs.dir directory.
        '''
        infiles = [i[0] for i in infiles]
        PipelinePeakcalling.mergeSortIndex(infiles, outfile)


elif PARAMS['IDR_poolinputs'] == "condition" and PARAMS['IDR_run'] != 1:
    @active_if(PARAMS['input'] != 0)
    @follows(mkdir('IDR_inputs.dir'))
    @split(filterInputBAMs, r'IDR_inputs.dir/*.bam')
    def makeIDRInputBams(infiles, outfiles):
        '''
        When IDR is going to be performed, inputs which are pooled by tissue
        and condition are automatically generated as these are always required.

        This function pools tissues and conditions when IDR is switched
        off if inputs pooled by condition are requested.

        The appropriate outputs from filterInputBAMs are identified and
        merged into a single BAM stored in the IDR_inputs.dir directory.
        '''
        outs = set(inputD.values())
        for out in outs:
            p = out.split("_")
            cond = p[0]
            tissue = p[1]

            # collect the appropriate bam files from their current location
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]
            innames = subdf['bamControl'].values
            innames = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in innames])
            out = "IDR_inputs.dir/%s" % out
            out = out.replace(".bam", "_filtered.bam")
            PipelinePeakcalling.mergeSortIndex(innames, out)


elif PARAMS['IDR_poolinputs'] == "condition" and PARAMS['IDR_run'] == 1:
    @active_if(PARAMS['input'] != 0)
    @follows(mkdir('IDR_inputs.dir'))
    @follows(mkdir('IDR_inputs.dir'))
    @transform(makePooledInputs, regex("IDR_inputs.dir/(.*).bam"),
               r'IDR_inputs.dir/\1.bam')
    def makeIDRInputBams(infiles, outfiles):
        '''
        If IDR is going to be run, pooled inputs are generated above so
        they don't need to be generated again if requested.
        '''
        pass


@follows(makeIDRInputBams)
@follows(filterInputBAMs)
@follows(makePooledBams)
@follows(makePooledInputs)
@follows(makePseudoBams)
@originate("peakcalling_bams_and_inputs.tsv")
def makeBamInputTable(outfile):
    '''
    Generates a tab delimited file - peakcalling_bams_and_inputs.tsv
    which links each filtered bam file in the peakcalling_bams.dir
    directory to the appropriate input in the IDR_inputs.dir
    directory.
    Uses the dictionary inputD generated as a global variable based
    on the user-specified design table plus pooled input files generated
    above.
    '''
    ks = inputD.keys()
    out = IOTools.openFile(outfile, "w")
    out.write('ChipBam\tInputBam\n')
    bamfiles = os.listdir("peakcalling_bams.dir")

    for k in ks:
        inputstem = inputD[k]
        chipstem = k
        chipstem = P.snip(chipstem)
        if PARAMS['input'] == 0:
            inputfile = "-"
        else:
            inputstem = P.snip(inputstem)
            inputfile = "IDR_inputs.dir/%s_filtered.bam" % inputstem

        for b in bamfiles:
            if b.startswith(chipstem) and b.endswith('bam'):
                out.write("peakcalling_bams.dir/%s\t%s\n" % (b, inputfile))
    out.close()


@transform(makeBamInputTable, suffix(".tsv"), ".load")
def loadBamInputTable(infile, outfile):
    P.load(infile, outfile)


@transform(makePseudoBams, suffix(".bam"), "_insertsize.tsv")
def estimateInsertSize(infile, outfile):
    '''
    Predicts insert size using MACS2 for single end data and using Bamtools
    for paired end data.
    Output is stored in insert_size.tsv
    '''
    PipelinePeakcalling.estimateInsertSize(infile, outfile,
                                           PARAMS['paired_end'],
                                           PARAMS['insert_alignments'],
                                           PARAMS['insert_macs2opts'])


@merge(estimateInsertSize, "insert_sizes.tsv")
def mergeInsertSizes(infiles, outfile):
    '''
    Combines insert size outputs into one file
    '''
    out = IOTools.openFile(outfile, "w")
    out.write("filename\tmode\tfragmentsize_mean\tfragmentsize_std\ttagsize\n")
    for infile in infiles:
        res = IOTools.openFile(infile).readlines()
        out.write("%s\t%s\n" % (infile, res[-1].strip()))
    out.close()


@transform(mergeInsertSizes, suffix(".tsv"), ".load")
def loadInsertSizes(infile, outfile):
    P.load(infile, outfile)


@follows(filtering)
@follows(loadInsertSizes)
@follows(loadDesignTable)
@follows(loadBamInputTable)
@transform(makePseudoBams, regex("(.*)_bams\.dir\/(.*)\.bam"),
           r"\1_bams.dir/\2.bam")
def preprocessing(infile, outfile):
    '''
    Dummy task to ensure all preprocessing has run and
    bam files are passed individually to the next stage.
    '''
    pass

#################################################################
# 3) Peakcalling
#################################################################


@follows(mkdir('macs2.dir'))
@transform(preprocessing,
           regex("peakcalling_bams.dir/(.*).bam"),
           add_inputs(makeBamInputTable),
           r"macs2.dir/\1.macs2")
def callMacs2peaks(infiles, outfile):
    '''
    Takes Bam and pairs with input using design files to
    call peaks using macs2

    Inputs
    ======
    bam file
    design file - looks up to identify which input file should be used
    for peakcalling
    instertsize.tsv - gets insert size to use for peak calling

    Output
    -----
    Macs2 output files
    hmmm- plus a couple of others - check the module file

    '''
    D = PipelinePeakcalling.readTable(infiles[1])
    bam = infiles[0]
    if PARAMS['input'] == 0:
        inputf = None
    else:
        inputf = D[bam]
    insertsizef = "%s_insertsize.tsv" % (P.snip(bam))

    peakcaller = PipelinePeakcalling.Macs2Peakcaller(
        threads=1,
        paired_end=PARAMS['paired_end'],
        tool_options=PARAMS['macs2_options'],
        tagsize=None)

    statement = peakcaller.build(bam, outfile,
                                 PARAMS['macs2_contigsfile'],
                                 inputf, insertsizef, PARAMS['IDR_run'],
                                 PARAMS['macs2_idrkeeppeaks'],
                                 PARAMS['macs2_idrsuffix'],
                                 PARAMS['macs2_idrcol'],
                                 PARAMS['macs2_broad_peak'])
    P.run()
    peakcaller.summarise(outfile)


@follows(mkdir('sicer_narrow.dir'))
@transform(preprocessing,
           regex("peakcalling_bams.dir/(.*).bam"),
           add_inputs(makeBamInputTable),
           r"sicer_narrow.dir/\1.narrow_sicer")
def callNarrowerPeaksWithSicer(infiles, outfile):
    '''
    Takes Bam and pairs with input using design files to
    call peaks using sicer

    Inputs
    ======
    bam file
    design file - looks up to identify which input file should be used
    for peakcalling
    instertsize.tsv - gets insert size to use for peak calling

    Output
    -----
    Sicer output files
    '''
    D = PipelinePeakcalling.readTable(infiles[1])
    bam = infiles[0]
    snip_bam = P.snip(bam)
    bam_name = snip_bam + "_insertsize.tsv"
    fragment_size = PipelinePeakcallingScrum.getFragment

    insert_size = DB.fetch_DataFrame("SELECT * FROM insert_sizes",
                                     PARAMS["database_name"])
    fragment_size = insert_size[insert_size['filename'].str.contains(bam_name)][
        'fragmentsize_mean']
    fragment_size = int(fragment_size.tolist()[0])

    window_size = PARAMS["sicer_narrow_window_size"]
    gap_size = PARAMS["sicer_narrow_gap_size"]
    fdr_threshold = PARAMS["sicer_fdr_threshold"]
    genome = PARAMS["genome"]
    redundancy_threshold = PARAMS["sicer_redundancy_threshold"]
    effective_genome_fraction = PARAMS['sicer_effective_genome_fraction']
    minfragsize = PARAMS['sicer_min_insert_size']
    maxfragsize = PARAMS['sicer_max_insert_size']

    # If there are no inputs
    if PARAMS['input'] == 0:
        inputf = None
    else:
        inputf = D[bam]

    peakcaller = PipelinePeakcalling.SicerPeakcaller(
        threads=1,
        tool_options=PARAMS['sicer_options'],
        window_size=window_size,
        gap_size=gap_size,
        fragment_size=fragment_size,
        fdr_threshold=fdr_threshold,
        effective_genome_fraction=effective_genome_fraction,
        genome=genome,
        redundancy_threshold=redundancy_threshold,
        minfragsize=minfragsize,
        maxfragsize=maxfragsize)

    statement = peakcaller.build(bam,
                                 outfile,
                                 controlfile=inputf,
                                 idr=PARAMS['IDR_run'],
                                 idrc=PARAMS['sicer_idrkeeppeaks'],
                                 idrcol=PARAMS['sicer_idrcol'],
                                 broad_peak=0)

    P.run()
    peakcaller.summarise(outfile, mode="narrow")


@follows(mkdir('sicer_broad.dir'))
@transform(preprocessing,
           regex("peakcalling_bams.dir/(.*).bam"),
           add_inputs(makeBamInputTable),
           r"sicer_broad.dir/\1.broad_sicer")
def callBroaderPeaksWithSicer(infiles, outfile):
    '''
    Takes Bam and pairs with input using design files to
    call peaks using sicer

    Inputs
    ======
    bam file
    design file - looks up to identify which input file should be used
    for peakcalling
    instertsize.tsv - gets insert size to use for peak calling

    Output
    -----
    Sicer output files
    '''
    D = PipelinePeakcalling.readTable(infiles[1])
    bam = infiles[0]
    snip_bam = P.snip(bam)
    bam_name = snip_bam + "_insertsize"
    insert_size = DB.fetch_DataFrame("SELECT * FROM insert_sizes",
                                     PARAMS["database_name"])
    fragment_size = insert_size[insert_size[
        'filename'].str.contains(bam_name)]['fragmentsize_mean']
    fragment_size = int(fragment_size.tolist()[0])

    window_size = PARAMS["sicer_broad_window_size"]
    gap_size = PARAMS["sicer_broad_gap_size"]
    fdr_threshold = PARAMS["sicer_fdr_threshold"]
    genome = PARAMS["genome"]
    redundancy_threshold = PARAMS["sicer_redundancy_threshold"]
    effective_genome_fraction = PARAMS['sicer_effective_genome_fraction']

    minfragsize = PARAMS['sicer_min_insert_size']
    maxfragsize = PARAMS['sicer_max_insert_size']

    # If there are no inputs
    if PARAMS['input'] == 0:
        inputf = None
    else:
        inputf = D[bam]

    peakcaller = PipelinePeakcalling.SicerPeakcaller(
        threads=1,
        tool_options=PARAMS['sicer_options'],
        window_size=window_size,
        gap_size=gap_size,
        fragment_size=fragment_size,
        fdr_threshold=fdr_threshold,
        effective_genome_fraction=effective_genome_fraction,
        genome=genome,
        redundancy_threshold=redundancy_threshold,
        minfragsize=minfragsize,
        maxfragsize=maxfragsize)

    statement = peakcaller.build(bam,
                                 outfile,
                                 controlfile=inputf,
                                 idr=PARAMS['IDR_run'],
                                 idrc=PARAMS['sicer_idrkeeppeaks'],
                                 idrcol=PARAMS['sicer_idrcol'],
                                 broad_peak=1)

    P.run()
    peakcaller.summarise(outfile, mode="broad")


@follows(callNarrowerPeaksWithSicer, callBroaderPeaksWithSicer)
def runSicer():
    pass

# list of peak callers to use
PEAKCALLERS = []
# list of peakcallers to use for IDR - currently IDR only works with a
# single peakcaller at a time
IDRPEAKCALLERS = []
# create dictionary of peakcallers and thier functions
mapToPeakCallers = {'macs2': (callMacs2peaks,),
                    'sicer': (runSicer,), }

# Call the peakcallers specified in the list
for x in P.asList(PARAMS['peakcalling_peakcallers']):
    PEAKCALLERS.extend(mapToPeakCallers[x])


@merge(PEAKCALLERS, "peakcalling_summary.tsv")
def summarisePeakCalling(infiles, outfile):
    bigtab = pd.DataFrame()
    for i in infiles:
        i = "%s_log.table" % i
        tab = pd.read_csv(i, sep="\t")
        bigtab = bigtab.append(tab)
    bigtab.to_csv(outfile, sep="\t", index=False)


@transform(summarisePeakCalling, suffix(".tsv"), ".load")
def loadPeakCallingStats(infile, outfile):
    P.load(infile, outfile)


@follows(loadPeakCallingStats)
def peakcalling():
    '''
    dummy task to collate upstream peakcalling tasks
    '''

################################################################
# 4) post peakcalling IDR Steps
################################################################

if PARAMS['IDR_run']:
    IDR_ON = True
else:
    IDR_ON = False


@active_if(IDR_ON)
@follows(peakcalling)
@follows(mkdir("peaks_for_IDR.dir"))
@transform(mapToPeakCallers[PARAMS['peakcalling_idrpeakcaller']],
           regex("(.*)/(.*)"),
           r"peaks_for_IDR.dir/\2.IDRpeaks")
def getIDRInputs(infile, outfile):
    '''
    Get the resulting peaks file from peakcalling
    and place them in IDR.dir so they can all be
    easilly found and indentified for IDR analysis

    inputs
    ======
    _IDRpeak files in peakcaller directorys
    (e.g. macs2.dir)

    output
    copy of _IDRpeak files in 'peaks_for_IDR.dir'
    '''
    IDRpeaks = "%s_IDRpeaks" % infile
    shutil.copy(IDRpeaks, outfile)


@active_if(IDR_ON)
@merge(getIDRInputs, "IDR_pairs.tsv")
def makeIDRPairs(infiles, outfile):
    '''
    generate table of files to pair up for
    IDR analysis

    inputs
    -----
    list of peak files in 'peaks_for_IDR.dir'

    Outputs
    -------
    table detailing the file pairings for IDR
    analysis
    '''
    useoracle = PARAMS['IDR_useoracle']
    PipelinePeakcalling.makePairsForIDR(infiles, outfile,
                                        PARAMS['IDR_useoracle'],
                                        df, submit=True)


@active_if(IDR_ON)
@transform(makeIDRPairs, suffix(".tsv"), ".load")
def loadIDRPairs(infile, outfile):
    P.load(infile, outfile)


@active_if(IDR_ON)
@follows(mkdir("IDR.dir"))
@split(makeIDRPairs, "IDR.dir/*.dummy")
def splitForIDR(infile, outfiles):
    '''
    infile = "IDR_pairs.tsv" file
    output =
    dummy file to act as placeholder for ruffus
    updated "IDR_pairs.tsv" file
    containainf tissue and condition information and
    the name of IDR output file
    '''
    pairs = pd.read_csv(infile, sep="\t")
    pairs['Condition'] = pairs['Condition'].astype('str')
    pairs['Tissue'] = pairs['Tissue'].astype('str')
    for p in pairs.index.values:
        p = pairs.ix[p]
        p1 = P.snip(p[0].split("/")[-1])
        p2 = P.snip(p[1].split("/")[-1])

        pairstring = "%s_v_%s" % (p1, p2)

        out = IOTools.openFile("IDR.dir/%s.dummy" % pairstring, "w")
        out.write("%s\n" % "\n".join(p))
        out.close()


@active_if(IDR_ON)
@transform(splitForIDR, suffix(".dummy"), ".tsv")
def runIDR(infile, outfile):
    ''' takes the  "IDR_pairs.tsv" detailing the files to be compared
    for IDR and uses this to run IDR analysis for the approriate files

    IDR_options = string from pipeline ini file detailing IDR options
    Different IDR comparisions (e.g. selfconistency, pooledconsistency or
    replicate consistancy might require different IDR thresholds) these can be
    set in the pipeline.ini file in the IDR section

    Oracle files = oracle peakset - see IDR analysis for details of what this
    means?

    '''
    lines = [line.strip() for line in IOTools.openFile(infile).readlines()]
    infile1, infile2, setting, oraclefile, condition, tissue = lines
    options = PARAMS['IDR_options']

    if setting == 'self_consistency':
        idrthresh = PARAMS['IDR_softthresh_selfconsistency']
        options += " %s" % PARAMS['IDR_options_selfconsistency']
    elif setting == "pooled_consistency":
        idrthresh = PARAMS['IDR_softthresh_pooledconsistency']
        options += " %s" % PARAMS['IDR_options_pooledconsistency']

    elif setting == "replicate_consistency":
        idrthresh = PARAMS['IDR_softthresh_replicateconsistency']
        options += " %s" % PARAMS['IDR_options_replicateconsistency']

    T = P.getTempFilename(".")
    statement = PipelinePeakcalling.buildIDRStatement(
        infile1, infile2,
        T,
        PARAMS['IDR_sourcecommand'],
        PARAMS['IDR_unsourcecommand'],
        idrthresh,
        idrPARAMS, options, oraclefile, test=True)

    P.run()
    lines = IOTools.openFile(T).readlines()
    os.remove(T)
    os.remove('%s.log' % T)

    if len(lines) >= 20:
        statement = PipelinePeakcalling.buildIDRStatement(
            infile1, infile2,
            outfile,
            PARAMS['IDR_sourcecommand'],
            PARAMS['IDR_unsourcecommand'],
            idrthresh,
            idrPARAMS, options, oraclefile)

        P.run()

    else:
        E.warn("""
        *******************************************************\
        IDR failed for %(infile1)s vs %(infile2)s - fewer than 20\
        peaks in the merged peak list\
        *******************************************************""" % locals())
        out = IOTools.openFile(outfile, "w")
        out.write("IDR FAILED - NOT ENOUGH PEAKS IN MERGED PEAK LIST")
        out.close()


@active_if(IDR_ON)
@transform(runIDR, suffix(".tsv"), ["_filtered.tsv",
                                    "_table.tsv"])
def filterIDR(infile, outfiles):
    '''
    Take the IDR output, which is in ENCODE narrowPeaks format if the input
    is narrowPeaks, gtf or bed and ENCODE broadPeaks format if the input is
    broadPeaks.
    Input is filtered based on whether it passes the soft IDR thresholds
    provided in the pipeline.ini.  Peaks which pass this threshold
    with have a score in the "globalIDR" column which is greater
    than -log(soft_threshold) where soft_threshold is the soft threshold
    provided in the pipeline.ini.
    Column headings are added and output is sorted by signalValue.
    '''
    IDRdata = pd.read_csv(infile, sep="\t", header=None)

    if 'FAILED' in IDRdata[0][0]:
        IDRdata.to_csv(outfiles[0], sep="\t")
        IDRpassed = 0
    else:
        IDRpassed = 1

        if idrPARAMS['idrsuffix'] == "broadPeak":
            IDRdata.columns = ["chrom", "chromStart", "chromEnd", "name",
                               "score", "strand", "signalValue",
                               "p-value", "q-value",
                               "localIDR", "globalIDR",
                               "rep1_chromStart", "rep2_chromEnd",
                               "rep1_signalValue", "rep2_chromStart",
                               "rep2_chromEnd", "rep2_signalValue"]
        else:
            IDRdata.columns = ["chrom", "chromStart", "chromEnd", "name",
                               "score", "strand", "signalValue", "p-value",
                               "q-value",
                               "summit", "localIDR", "globalIDR",
                               "rep1_chromStart", "rep2_chromEnd",
                               "rep1_signalValue", "rep1_summit",
                               "rep2_chromStart", "rep2_chromEnd",
                               "rep2_signalValue", "rep2_summit"]

        IDRdataP = IDRdata[IDRdata['score'] == 1000]
        IDRdataF = IDRdata[IDRdata['score'] != 1000]

        IDRdataP = IDRdataP.sort_values('signalValue', ascending=False)
        IDRdataF = IDRdataF.sort_values('signalValue', ascending=False)
        IDRdataP.to_csv(outfiles[0], sep="\t")

    H = ['Total_Peaks', 'Peaks_Passing_IDR', 'Peaks_Failing_IDR',
         'Percentage_Peaks_Passing_IDR', 'IDR_Successful']

    if IDRpassed == 1:
        T = ((len(IDRdata), len(IDRdataP), len(IDRdataF),
              round(float(len(IDRdataP)) / float(len(IDRdata)), 4) * 100,
              "TRUE"))
    else:
        T = ((0, 0, 0, 0, 0, "FALSE"))

    out = IOTools.openFile(outfiles[1], "w")
    out.write("%s\n" % "\t".join(H))
    out.write("%s\n" % "\t".join([str(t) for t in T]))


@active_if(IDR_ON)
@merge((filterIDR, makeIDRPairs), "IDR_results.tsv")
def summariseIDR(infiles, outfile):
    '''
    '''
    pooledc, selfc, repc = (PARAMS['IDR_softthresh_pooledconsistency'],
                            PARAMS['IDR_softthresh_selfconsistency'],
                            PARAMS['IDR_softthresh_replicateconsistency'])

    PipelinePeakcalling.summariseIDR(infiles, outfile, pooledc, selfc, repc)


@active_if(IDR_ON)
@transform(summariseIDR, suffix(".tsv"), ".load")
def loadIDRsummary(infile, outfile):
    P.load(infile, outfile)


@active_if(IDR_ON)
@transform(summariseIDR, suffix("results.tsv"), "QC.tsv")
def runIDRQC(infile, outfile):
    PipelinePeakcalling.doIDRQC(infile, outfile)


@active_if(IDR_ON)
@transform(runIDRQC, suffix(".tsv"), ".load")
def loadIDRQC(infile, outfile):
    P.load(infile, outfile)


@active_if(IDR_ON)
@follows(mkdir("conservative_peaks.dir"))
@split(summariseIDR, "conservative_peaks.dir\/*\.tsv")
def findConservativePeaks(infile, outfiles):
    '''function selects row from IDR_results.tsv that represents the
    conservative peak list'''
    tab = pd.read_csv(infile, sep="\t")
    cps = tab[tab['Conservative_Peak_List'] == 'True']
    experiments = cps['Experiment'].values
    peakfiles = cps['Output_Filename'].values

    peakfiles = ["IDR.dir/%s" % i.replace("_table", "") for i in peakfiles]
    i = 0
    for peakfile in peakfiles:
        outnam = "conservative_peaks.dir/%s.tsv" % experiments[i]
        PipelinePeakcalling.makeLink(peakfile, outnam)
        i += 1


@active_if(IDR_ON)
@follows(mkdir("optimal_peaks.dir"))
@split(summariseIDR, "conservative_peaks.dir\/*\.tsv")
def findOptimalPeaks(infile, outfiles):
    '''function selects row from IDR_results.tsv that represents the
    optimal peak list'''
    tab = pd.read_csv(infile, sep="\t")
    cps = tab[tab['Optimal_Peak_List'] == 'True']
    experiments = cps['Experiment'].values
    peakfiles = cps['Output_Filename'].values

    peakfiles = ["IDR.dir/%s" % i.replace("_table", "") for i in peakfiles]

    i = 0
    for peakfile in peakfiles:
        outnam = "optimal_peaks.dir/%s.tsv" % experiments[i]
        PipelinePeakcalling.makeLink(peakfile, outnam)
        i += 1


@active_if(IDR_ON)
@follows(loadIDRPairs)
@follows(loadIDRsummary)
@follows(loadIDRQC)
@follows(findConservativePeaks)
@follows(findOptimalPeaks)
@follows(runIDRQC)
def IDR():
    pass


################################################################
# QC Steps


@merge(("design.tsv", makeBamInputTable),
       ["ChIPQC_design_conservative.tsv",
       "ChIPQC_design_optimal.tsv"])
def makeCHIPQCInputTables(infiles, outfiles):
    design = pd.read_csv(infiles[0], sep="\t")
    inputs = pd.read_csv(infiles[1], sep="\t")
    inputs['SampleID'] = inputs[
        'ChipBam'].str.split("/").str.get(-1).str.split(".").str.get(0)
    inputs['SampleID'] = [i.replace("_filtered", "")
                          for i in inputs['SampleID']]
    tab = design.merge(inputs)
    tab = tab.drop("ControlID", 1)
    tab = tab.drop("bamReads", 1)
    tab = tab.rename(columns={"ChipBam": "bamReads"})
    tab = tab[['SampleID', 'Tissue', 'Condition', 'Replicate',
               'bamReads']]
    tab = tab.rename(columns={'Condition': 'Factor'})
    tab['Factor'] = tab['Factor'].astype('str')
    tab['Tissue'] = tab['Tissue'].astype('str')

    tab['Peaks'] = ("conservative_peaks.dir/" + tab['Factor'] + "_" +
                    tab['Tissue'] + ".tsv")

    tab.to_csv(outfiles[0], sep="\t", index=None)

    tab['Peaks'] = ("optimal_peaks.dir/" + tab['Factor'] + "_" +
                    tab['Tissue'] + ".tsv")
    tab.to_csv(outfiles[1], sep="\t", index=None)

# TODO
# @follows(mkdir("ChIPQC.dir"))
# @transform(makeCHIPQCInputTable,regex("(.*)_(.*).tsv"), r'ChIPQC.dir/\1.pdf')
# def runCHIPQC(infiles, outfiles):
#    R('''''')


@follows(filtering, peakcalling, IDR)
def full():
    ''' runs entire pipeline '''
    pass


###############################################################
# Report functions
###############################################################


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


@follows(mkdir("%s/bamfiles" % PARAMS["web_dir"]),
         mkdir("%s/medips" % PARAMS["web_dir"]),
         )
def publish():
    '''publish files to web directory'''

    # directory : files

    # publish web pages
    # P.publish_report(export_files=export_files)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
