"""
pipeline_windows.py - Window based genomic analysis
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline takes mapped reads from ChIP-Seq experiments
such has chromatin marks, MeDIP and performs analyses
of the genomic read distribution. This contrasts with
:doc:`pipeline_intervals`, which annotates a set of
non-overlapping intervals.

The pipeline performs the following analyses:

Window based analysis
    The pipeline defines windows across the genome
    add counts the reads mapping into the windows.
    It then detects if there are any differences
    in window read counts between different experimental
    conditions

Genomic context analysis
    The genome is divided into annotations such
    as repeat, exon, ..... It then detects if there
    are any differences in these regions.

Meta-gene profiling
    Compute read distributions across genes.

Enrichment analysis
    For genomic context, the pipeline detects
    if tags lie more frequently then expected inside
    particular genomic contexts.



Methods
=======

Window based analysis
---------------------

   1. Identify differentially occupied regions
   4. Filter DMRs
   5. Calculate DMR statistics
   6. Produce report (SphinxReport)


Tiling strategies
-----------------

The pipeline implements different tiling strategies.

variable width
   variable width tiles. Tiles are defined based on regions that contain
   short reads and are present in a minimum number of samples.

fixwidth_nooverlap
   tiles of size ``tiling_window_size`` with adjacent tiles not overlapping.

fixwidth_overlap
   tiles of size ``tiling_window_size`` with adjacent tiles overlapping by
   by 50%.

cpg
   windows of size ``tiling_window_size`` are defined around CpG sites.
   Overlapping windows are merged and only windows with a minimum number
   (``tiling_min_cpg``) of CpG sites are kept.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----
The input is mapped reads in bam files in the working directory or
linked from the working directory.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.bam


.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

Pipeline output
===============

Requirements:

* bedtools >= 2.21.0
* ucsctools

Code
====

"""

# load modules
from ruffus import transform, merge, mkdir, follows, \
    regex, suffix, add_inputs, collate

import sys
import os
import re
import glob
import csv
import numpy
import sqlite3
import pandas

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineWindows as PipelineWindows
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMappingQC as PipelineMappingQC

from rpy2.robjects import r as R

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    prefix="annotations_",
    update_interface=True))


###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.AutoSample

METHODS = P.asList(PARAMS["methods"])


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


# @P.add_doc(PipelineWindows.convertReadsToIntervals)
@follows(mkdir("tags.dir"))
@transform('*.bam',
           regex("(.*).bam"),
           r"tags.dir/\1.bed.gz")
def prepareTags(infile, outfile):
    '''
    Parameters
    ----------
    infile: str
        Filename of input file in :term:`bam` format
    outfile: str
        Filename of output file in :term:`bed` format.
    filtering_quality : int
        :term:`PARAMS`
        If set, remove reads with a quality score below given threshold.
    filtering_dedup : bool
        :term:`PARAMS`
        If True, deduplicate data.
    filtering_dedup_method : string
        :term:`PARAMS`
        Deduplication method. Possible options are ``picard`` and
        ``samtools``.
    filtering_nonunique : bool
        :term:`PARAMS`
        If True, remove non-uniquely matching reads.
    '''
    PipelineWindows.convertReadsToIntervals(
        infile,
        outfile,
        filtering_quality=PARAMS.get('filtering_quality', None),
        filtering_dedup='filtering_dedup' in PARAMS,
        filtering_dedup_method=PARAMS['filtering_dedup_method'],
        filtering_nonunique=PARAMS.get('filtering_nonunique', False))


# @P.add_doc(PipelineWindows.countTags)
@transform(prepareTags, suffix(".bed.gz"), ".tsv")
def countTags(infile, outfile):
    PipelineWindows.countTags(infile, outfile)


@merge(countTags, "tag_counts.load")
def loadTagCounts(infiles, outfile):
    '''Load tag counts representing read counts from bed files into database
       as table tag_counts

       Parameters
       ----------
       infiles: list
           filenames of :term:`tsv` formatted files containing tag counts
       outfile: str
           filename of database loading logfile.
       '''
    P.mergeAndLoad(infiles, outfile, columns=(0, 2),
                   suffix=".tsv")


# @P.add_doc(PipelineMappingQC.loadPicardDuplicateStats)
@merge(prepareTags, "picard_duplicates.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite
    table picard_duplicates.
    '''
    PipelineMappingQC.loadPicardDuplicateStats(
        infiles, outfile, pipeline_suffix=".bed.gz")


@follows(mkdir("background.dir"))
@transform("*[Ii]nput*.bw",
           regex("(.*).bw"),
           r"background.dir/\1.bed.gz")
def buildBackgroundWindows(infile, outfile):
    '''compute regions with high background count in each input (untreated)
    bigwig file.  Bigwigs can be generated with pipeline_mapping or
    bam2bigwig.py.

    Parameters
    ----------
    infile: str
        filename of :term:`bigwig` file showing coverage in input

    genome_dir: str
        :term:`PARAMS`
        path to indexed genome

    filtering_background_density: int
        :term:`PARAMS`
        regions above this threshold are classed as high background counts
        and are removed.

    outfile: str
        filename of output :term:`bed` file showing high coverage regions
    '''

    job_memory = "16G"

    statement = '''
    cgat wig2bed
             --bigwig-file=%(infile)s
             --genome-file=%(genome_dir)s/%(genome)s
             --threshold=%(filtering_background_density)f
             --method=threshold
             --log=%(outfile)s.log
    | bgzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


@merge(buildBackgroundWindows, "background.dir/background.bed.gz")
def mergeBackgroundWindows(infiles, outfile):
    '''Build a single bed file of regions with elevated background by
    merging different input samples / runs.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`bed` formatted files containing
        regions of elevated background.

    filtering_background_extension: int
        :term:`PARAMS`
        number of bases to extend elevated regions from read mapping positions.

    outfile: str:
        filename for combined :term:`bed` file
    '''

    if len(infiles) == 0:
        # write a dummy file with a dummy chromosome
        # an empty background file would otherwise cause
        # errors downstream in bedtools intersect
        outf = IOTools.openFile(outfile, "w")
        outf.write("chrXXXX\t1\t2\n")
        outf.close()
        return

    infiles = " ".join(infiles)
    genomefile = os.path.join(PARAMS["annotations_dir"],
                              PARAMS["annotations_interface_contigs"])
    statement = '''
    zcat %(infiles)s
    | bedtools slop -i stdin
                -b %(filtering_background_extension)i
                -g %(genomefile)s
    | sort -k 1,1 -k2,2n
    | bedtools merge -i -
    | bgzip
    > %(outfile)s
    '''

    P.run()


@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS["annotations_interface_cpg_bed"]),
           regex(".*/([^/]*).bed.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_genomic_context_bed"])),
           "cpg_context.tsv.gz")
def buildCpGAnnotation(infiles, outfile):
    '''annotate the location of CpGs using a bed file showing
    CpG regions and a bed file showing genomic context - counts the number of
    CpGs overlapping different genome annotations.

    Parameters
    ----------
    infiles: list
       consists of
    infiles[0]: str
       :term:`bed` annotation file showing locations of CpGs
    infiles[1]: str
       :term:`bed` annotation file showing genomic context
    outfile: str
       :term:`tsv` formatted file showing how many CpGs overlap each genomic
        annotation type.
    '''
    cpg_bed, context_bed = infiles
    statement = '''
    cgat bam_vs_bed
           --min-overlap=0.5 %(cpg_bed)s %(context_bed)s
    | gzip
    > %(outfile)s'''

    P.run()


@transform(buildCpGAnnotation, suffix(".tsv.gz"), ".load")
def loadCpGAnnotation(infile, outfile):
    '''Load CpG genomic context data into database
       as table cpg_context

       Parameters
       ----------
       infile: str
           filename of :term:`tsv` file containing CpG genomic context
       outfile: str
           filename of database loading logfile.
       '''
    P.load(infile, outfile)


@transform(prepareTags, suffix(".bed.gz"), ".covered.bed.gz")
def buildCoverageBed(infile, outfile):
    '''
    Build a :term:`bed` file of regions covered by reads.
    Intervals containing only few tags (tiling_min_reads) are removed.

    Parameters
    ----------
    infile: str
        filename of :term:`bed` file containing tag counts

    medips_extension: int
        :term:`PARAMS`
        reads are extended to represent the fragment size from the sonication
        step of MEDIP analysis by this number of bases.

    tiling_min_reads: int
        :term:`PARAMS`
        minimum number of reads to class a region as covered by reads

    outfile: str
        filename of output :term:`bed` file to show regions covered by more
        than tiling_min_reads reads.
    '''

    statement = '''
    zcat %(infile)s
    | cut -f 1,2,3
    | cgat bed2bed
          --method=merge
          --merge-distance=%(medips_extension)i
          --log=%(outfile)s.log
          --merge-min-intervals=%(tiling_min_reads)i
    | gzip
    > %(outfile)s
    '''
    P.run()


@transform(buildCoverageBed, suffix(".bed.gz"), ".tsv.gz")
def buildCpGComposition(infile, outfile):
    '''
    Compute CpG density across regions covered by reads

    Parameters
    ----------
    infile: str
        filename of :terms:`bed` file showing regions covered by reads

    genome_dir: str
        :term:`PARAMS`
        directory containing indexed reference genome

    genome: str
        :term:`PARAMS`
        name of reference genome file

    outfile: str
        filename of :terms:`tsv` file to write the CpG density information
    '''

    statement = '''
    zcat %(infile)s
    | cgat bed2table
    --counter=composition-cpg
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s
    '''
    P.run()


@merge(buildCoverageBed, "tags.dir/genomic.covered.tsv.gz")
def buildReferenceCpGComposition(infiles, outfile):
    '''
    Compute CpG densities across reference windows across
    the genome.

    This will take the first file of the input and
    shuffle the intervals, and then compute.

    Using fixed size windows across the genome results in
    a very discretized distribution compared to the other
    read coverage tracks which have intervals of different size.

    Parameters
    ----------
    infiles: list
        list of filenames of bed files of regions covered by reads

    annotations_dir: str
        :term:`PARAMS`
        directory containing annotations info

    annotations_interface_contigs: str
        :term:`PARAMS`
        filename of contig size annotation

    annotations_interface_gaps_bed: str
        :term:`PARAMS`
        filename of annotation of gaps

    outfile: str
        filename to write the reference CpG composition in :term:`tsv` format
    '''

    infile = infiles[0]
    contig_sizes = os.path.join(PARAMS["annotations_dir"],
                                PARAMS["annotations_interface_contigs"])
    gaps_bed = os.path.join(PARAMS["annotations_dir"],
                            PARAMS["annotations_interface_gaps_bed"])

    # remove windows which are more than 50% N - column 17
    statement = '''bedtools shuffle
                      -i %(infile)s
                      -g %(contig_sizes)s
                      -excl %(gaps_bed)s
                      -chromFirst
    | cgat bed2table
          --counter=composition-cpg
          --genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s
    '''
    P.run()

    # | awk '$1 !~ /%(tiling_remove_contigs)s/'
    # | awk '$1 == "contig" || $17 < 0.5'


@transform((buildCpGComposition,
            buildReferenceCpGComposition),
           suffix(".tsv.gz"), ".cpghist.tsv.gz")
def histogramCpGComposition(infile, outfile):
    '''Build histogram of CpG density in all regions covered by reads.

    Parameters
    ----------
    infile: tuple
       filenames of CpG composition of regions covered by reads (infile[0])
       and permuted reference version (infile[1])
    outfile:
       filename for histogram in :term:`tsv` format
    '''

    statement = '''
    zcat %(infile)s
    | cgat csv_cut pCpG
    | cgat data2histogram --bin-size=0.01
    | gzip
    > %(outfile)s
    '''
    P.run()


@merge(histogramCpGComposition, "pcpg_in_coveredregions.load")
def loadCpgCompositionHistogram(infiles, outfile):
    '''Load histograms of CpG Density in regions covered by reads into
    database table - pcpg_in_coveredregions

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`tsv` formatted histograms showing CpG
        composition
    outfile: str
        filename for database load logfile
    '''

    P.mergeAndLoad(infiles, outfile,
                   regex="/(.*).cpghist.tsv.gz",
                   row_wise=False)


@transform((buildCpGComposition, buildReferenceCpGComposition),
           suffix(".tsv.gz"),
           ".composition.load")
def loadCpGComposition(infile, outfile):
    '''Load CpG genomic composition data into database
       as table genomic.covered.composition

       Parameters
       ----------
       infile: tuple
          filenames of CpG composition of regions covered by reads (infile[0])
          and permuted reference version (infile[1])
       outfile: str
           filename of database loading logfile.
       '''
    P.load(infile, outfile)


@transform(prepareTags,
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS["annotations_interface_cpg_bed"])),
           ".cpg_coverage.gz")
def buildCpGCoverage(infiles, outfile):
    '''Count the number of times each CpG is covered by reads.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        filename of bed file containing read counts
    infiles[1]: str
        filename of :term:`bed` file containing CpG annotation
    outfile: str
        filename for histogram showing how many times CpGs are covered by reads

    ??Reads are processed in the same way as by buildCoverageBed.

    '''

    # coverageBed is inefficient. If bedfile and cpgfile
    # were sorted correspondingly the overlap analysis
    # could be done in very little memory.

    infile, cpg_file = infiles

    job_memory = "32G"

    statement = '''
    zcat %(infile)s
    | bedtools coverage -a stdin -b %(cpg_file)s -counts
    | cut -f 6
    | cgat data2histogram
    | gzip
    > %(outfile)s
     '''
    P.run()


@merge(buildCpGCoverage, "cpg_coverage_by_reads.load")
def loadCpGCoverage(infiles, outfile):
    '''load cpg coverage data - number of reads covering a CpG.'''
    P.mergeAndLoad(infiles, outfile,
                   regex="/(.*).cpg_coverage.gz",
                   row_wise=False)


@follows(loadCpGCoverage, loadCpGComposition, loadCpGAnnotation)
def gc():
    pass


@merge((buildCoverageBed, mergeBackgroundWindows), "windows.bed.gz")
def buildWindows(infiles, outfile):
    '''Build tiling windows according to parameter tiling_method.  Windows in
    background are removed.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        filename of :term:`bed` file showing regions covered by reads
    infiles[1]: str
        filename of :term:`bed` file showing regions in background (to remove)
    tiling_method: str
        :terms:`PARAMS`
        can be fixwidth_overlap, fixwidth_nooverlap, varwidth, cpg or use
        a bed file.
    tiling_window_size: int
        :terms:`PARAMS`
        window size for fixed width windows
    tiling_min_cpg: int
        minimum number of CpGs for CpG tiling method
    tiling_remove_contigs: str
        patterns to match for contigs to remove
    outfile: str
        filename of :terms:`bed` file to show the windows

    '''

    tiling_method = PARAMS["tiling_method"]

    coverage_bed, background_bed = infiles[:-1], infiles[-1]

    coverage_bed = " ".join(coverage_bed)

    if tiling_method == "varwidth":

        infiles = " ".join(infiles)

        statement = '''
        zcat %(coverage_bed)s
        | sort -k1,1 -k2,2n
        | cgat bed2bed
              --method=merge
              --merge-distance=0
              --log=%(outfile)s.log
        '''

    elif tiling_method == "fixwidth_nooverlap":

        statement = '''cgat genome_bed
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_window_size)i
                      --shift-size=%(tiling_window_size)i
                      --log=%(outfile)s.log'''

    elif tiling_method == "fixwidth_overlap":

        assert PARAMS["tiling_window_size"] % 2 == 0
        shift = PARAMS["tiling_window_size"] // 2

        statement = '''cgat genome_bed
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_window_size)i
                      --shift-size=%(shift)i
                      --log=%(outfile)s.log'''

    elif tiling_method == "cpg":

        statement = '''cat %(genome_dir)s/%(genome)s.fasta
                       | cgat fasta2bed
                      --method=windows-cpg
                      --window-size=%(tiling_window_size)i
                      --min-cpg=%(tiling_min_cpg)i
                      --log=%(outfile)s.log'''

    elif os.path.exists(tiling_method):
        # existing file
        statement = '''mergeBed -i %(tiling_method)s'''

    else:
        raise ValueError("unknow tiling method '%s'" % tiling_method)

    statement += '''
        | awk '$1 !~ /%(tiling_remove_contigs)s/'
        | bedtools intersect -v -wa -a stdin -b %(background_bed)s
        | gzip
        > %(outfile)s
    '''

    P.run()


@transform(buildWindows,
           suffix(".bed.gz"),
           ".stats")
def buildWindowStats(infile, outfile):
    '''Compute tiling window size statistics from bed file.

    Parameters
    ----------
    infile: str
        filename of :term:`bed` file representing tiled genome

    outfile: str
        filename of :term:`tsv` file for window size statistics as histograms
    '''

    statement = '''
    zcat %(infile)s
    | cgat gff2histogram
                   --force-output
                   --format=bed
                   --output-section=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()


@transform(buildWindowStats,
           suffix(".stats"),
           "_stats.load")
def loadWindowStats(infile, outfile):
    '''
    Load window size histograms to database table - <track>_stats, where
    track is the prefix of the input bam filename.

    Parameters
    ----------
    infile: str
        name of :term:`csv` file with window size histograms
    outfile: str
        filename of database load log file.
    '''
    P.load(infile + ".hist.tsv", P.snip(infile, ".stats") + "_hist" + ".load")
    P.load(infile + ".stats.tsv", outfile)


@transform(buildWindows,
           suffix(".bed.gz"),
           ".composition.tsv.gz")
def buildWindowComposition(infile, outfile):
    '''
    Compute length, cpg composition and genomic context of windows

    Parameters
    ----------
    infile: str
        filename of :term:`bed` file containing window positions
    outfile: str
        filename of :term:`tsv` formatted file to write statistics.
    '''
    statement = '''
    zcat %(infile)s
    | cgat bed2table
    --log=%(outfile)s.log
    --counter=length
    --counter=composition-cpg
    --counter=composition-na
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s
    '''
    P.run()


@transform(buildWindows,
           suffix(".bed.gz"),
           ".bigbed")
def buildBigBed(infile, outfile):
    '''
    Builds a :term:`bigBed` file with intervals that are covered by reads in
    each experiment.

    Parameters
    ----------
    infile: str
        filename of :term:`bed` file showing which intervals are covered by
        reads

    annotations_interface_contigs: str
        filename of annotation of contig lengths

    outfile: str
        filename of :term:`bigbed` file to output the results
    '''

    tmpfile = P.getTempFilename()

    contig_sizes = os.path.join(
        PARAMS["annotations_dir"], PARAMS["annotations_interface_contigs"])

    statement = '''
    zcat %(infile)s > %(tmpfile)s;
    bedToBigBed %(tmpfile)s %(contig_sizes)s %(outfile)s;
    rm -f %(tmpfile)s
    '''
    P.run()

    try:
        os.unlink(tmpfile)
    except OSError:
        pass


# @P.add_doc(PipelineWindows.countTagsWithinWindows
@follows(mkdir("counts.dir"))
@transform(prepareTags,
           regex(".*/(.*).bed.gz"),
           add_inputs(buildWindows),
           r"counts.dir/\1.counts.bed.gz")
def countTagsWithinWindows(infiles, outfile):
    '''
    Count the number of reads mapped to each window

    Parameters
    ----------
    infiles: list
        list of filenames

    infiles[0]: str
        filename of :term:`bed` format file containing read counts for
        all sites

    infiles[1]: str
        filename of :term: `bed` format file containing window positions

    tiling_counting_method: str
        :term:`PARAMS`
        can be "midpoint" or "nucleotide"
        midpoint counts the number of reads overlapping the midpoint of the
        window by at least one base
        nucleotide counts the number of reads overlapping the window by at
        least one base.

    outfile: str
        filename for :term:`bed` formatted file of read counts per window
    '''
    bedfile, windowfile = infiles
    PipelineWindows.countTagsWithinWindows(
        bedfile,
        windowfile,
        outfile,
        counting_method=PARAMS['tiling_counting_method'],
        job_memory=PARAMS['tiling_counting_memory'])


# @P.add_doc(PipelineWindows.aggregateWindowsTagCounts)
@merge(countTagsWithinWindows,
       r"counts.dir/windows_counts.tsv.gz")
def aggregateWindowsTagCounts(infiles, outfile):
    '''
    Aggregate tag counts into a single file.

    Parameters
    ----------
    infiles: list
        filenames of all :term:`bed` formatted window read count files

    outfile: str
        output filename for compiled window read counts
    '''

    PipelineWindows.aggregateWindowsTagCounts(infiles,
                                              outfile,
                                              regex="(.*).counts.bed.gz")


# @P.add_doc(PipelineWindows.countTagsWithinWindows)
@follows(mkdir('contextstats.dir'))
@transform(prepareTags,
           regex(".*/(.*).bed.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_genomic_context_bed"])),
           r"contextstats.dir/\1.counts.bed.gz")
def countTagsWithinContext(infiles, outfile):
    '''collect context stats of BED files.

    Examines the genomic context to where tags are located.

    A tag is assigned to the genomic context that it overlaps by at
    least 50%. Long tags that map across several non-overlapping
    contexts might be dropped.

    Parameters
    ----------
    infiles: str
        list of filenames

    infiles[0]: str
        filename of :term:`bed` file containing tag counts

    infiles[1]: str
        filename of :term:`bed` file containing genomic context

    outfile: str
        filename for :term:`bed` formatted file of genomic context of tags

    '''
    tagfile, windowfile = infiles
    PipelineWindows.countTagsWithinWindows(tagfile,
                                           windowfile,
                                           outfile,
                                           counting_method="midpoint",
                                           job_memory="4G")


# @P.add_doc(PipelineWindows.aggregateWindowsTagCounts)
@merge(countTagsWithinContext,
       r"contextstats.dir/context_counts.tsv.gz")
def aggregateContextTagCounts(infiles, outfile):
    '''Aggregate tag counts into a single file.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`bed` files containing tag counts for
        genomic contexts
    outfile: str
        filename for :term:`tsv` formatted file showing tag counts for
        genomic contexts for all samples.

    '''
    PipelineWindows.aggregateWindowsTagCounts(infiles,
                                              outfile,
                                              regex="(.*).counts.bed.gz")


@transform((aggregateWindowsTagCounts, aggregateContextTagCounts),
           suffix(".tsv.gz"),
           "_normed.tsv.gz")
def normalizeTagCounts(infile, outfile):
    '''
    Normalises tag counts from different experiments to make them comparable.

    Parameters
    ----------
    infile: str
        filename of file containing tag counts across windwos for all
        samples in :term:`tsv` format
    tags_normalization_method: str
        can be deseq-size factors, total-column, total-row, total-count
        deseq-size-factors - use normalisation implemented in DEseq
        total-column - divide counts by column total
        total-row - divide counts by the value in a row called 'total'
        total-count - normalised all values in column by the ratio of the
        per column sum of counts and the average column count
        across all rows.
    outfile: str
        :term:`tsv` filename to write the normalised tag counts.

    '''
    PipelineWindows.normalizeTagCounts(
        infile,
        outfile,
        method=PARAMS["tags_normalization_method"])


@transform((aggregateWindowsTagCounts, aggregateContextTagCounts),
           suffix(".tsv.gz"), ".load")
def loadWindowsTagCounts(infile, outfile):
    '''
    Load a sample of window composition data and context data
    into a database for QC purposes - generates two tables,
    context_counts and windows_counts.

    Parameters
    ----------
    infile: str
        filename of aggregated window or context tag counts
    outfile: str
        logfile for database load
    '''
    P.load(infile, outfile, limit=10000, shuffle=True)


def getInput(track):
    '''return a list of input tracks associated with track.

    Associations can be defined in the .ini file in the section
    [input]. For example, the following snippet associates track
    track1 with the bamfiles :file:`track1.bam` and :file:`track2.bam`::

       [input]
       track1.bam=input1.bam,input2.bam

    Glob expressions are permitted.

    Default tracks can be specified using a placeholder ``%``. The
    following will associate all tracks with the same bam file::

        [bams]
        %=all.bam

    Parameters
    ----------
    track: str
        filename of bam file of interest

    '''

    input_files = []

    # configparser by default converts option names to lower case
    fn = track.asFile()
    fn = fn.lower()

    if "input_%s" % fn in PARAMS:
        input_files.extend(P.asList(PARAMS["input_%s" % fn]))
    elif P.CONFIG.has_section("input"):
        for pattern, value in P.CONFIG.items("input"):
            if "%" in pattern:
                pattern = re.sub("%", "\S+", pattern)
            if re.search(pattern, fn):
                input_files.extend(P.asList(value))
    input_files = [re.sub('[^0-9a-zA-Z]+', '_', i) for i in input_files]
    return input_files


def mapTrack2Input(tracks):
    '''
    Given a list of tracks, return a dictionary mapping a track to its input

    Parameters
    ----------
    tracks: list
        list of strings representing filenames of :term:`bam` file
    '''

    # select columns in foreground and background
    map_track2input = {}
    for idx, track in enumerate(tracks):

        if track == "interval_id":
            continue

        try:
            t = Sample(tablename=track)
        except ValueError as msg:
            print(msg)
            continue

        input_files = getInput(t)

        # currently only implement one input file per track
        assert len(input_files) <= 1, "%s more than input: %s" % (
            track, input_files)

        if len(input_files) == 0:
            map_track2input[track] = None
        else:
            map_track2input[track] = Sample(filename=input_files[0]).asTable()

    return map_track2input


@transform(loadWindowsTagCounts,
           suffix(".load"),
           "_l2foldchange_input.tsv.gz")
def buildWindowsFoldChangesPerInput(infile, outfile):
    '''
    Compute fold changes for each sample compared to appropriate input.
    If no input is present, simply divide by average.

    Parameters
    ----------
    infile: str
        filename of log from database load of windows tag counts
    outfile: str
        filename of :term:`tsv` formatted file to write the fold change data

    '''

    # get all data
    dbh = connect()
    cc = dbh.cursor()
    cc.execute("SELECT * FROM windows_counts")
    data = cc.fetchall()

    # transpose, remove interval_id column
    data = list(zip(*data))
    columns = [x[0] for x in cc.description]

    map_track2input = mapTrack2Input(columns)
    take_tracks = [x for x, y in enumerate(columns) if y in map_track2input]
    take_input = [x for x, y in enumerate(
        columns) if y in list(map_track2input.values()) and y is not None]

    # build data frame
    dataframe = pandas.DataFrame(
        dict([(columns[x], data[x]) for x in take_tracks]))
    dataframe = dataframe.astype('float64')
    dataframe_input = pandas.DataFrame(
        dict([(columns[x], data[x]) for x in take_input]))

    # add pseudocounts
    pseudocount = 1
    for column in dataframe.columns:
        dataframe[column] += pseudocount
    for column in dataframe_input.columns:
        dataframe_input[column] += pseudocount

    # compute normalization ratios
    # total_input / total_column
    ratios = {}
    for column in dataframe.columns:
        i = map_track2input[column]
        if i is not None:
            ratios[column] = dataframe_input[
                i].median() / dataframe[column].median()
        else:
            ratios[column] = None

    for column in dataframe.columns:
        if ratios[column] is not None:
            # normalize by input
            dataframe[column] *= ratios[column] / \
                dataframe_input[map_track2input[column]]
        else:
            # normalize by median
            dataframe[column] /= dataframe[column].median()

    dataframe = numpy.log2(dataframe)

    dataframe.to_csv(IOTools.openFile(outfile, "w"),
                     sep="\t", index=False)


@transform(loadWindowsTagCounts,
           suffix(".load"),
           "_l2foldchange_median.tsv.gz")
def buildWindowsFoldChangesPerMedian(infile, outfile):
    '''
    Compute l2fold changes for each sample compared to the median count
    in the sample.

    Parameters
    ----------
    infile: str
        name of database load logfile for :term:`tsv` file of read counts
        across windows
    outfile: str
        filename for :term:`tsv` formatted file to write the fold change data
    '''

    # get all data
    dbh = connect()
    cc = dbh.cursor()
    cc.execute("SELECT * FROM windows_counts")
    data = cc.fetchall()

    # transpose, remove interval_id column
    data = list(zip(*data))
    columns = [x[0] for x in cc.description]

    take_tracks = [x for x, y in enumerate(columns) if y != "interval_id"]
    # build data frame
    dataframe = pandas.DataFrame(
        dict([(columns[x], data[x]) for x in take_tracks]))
    dataframe = dataframe.astype('float64')

    # add pseudocounts
    pseudocount = 1
    for column in dataframe.columns:
        dataframe[column] += pseudocount

    for column in dataframe.columns:
        dataframe[column] /= dataframe[column].median()

    dataframe = numpy.log2(dataframe)

    dataframe.to_csv(IOTools.openFile(outfile, "w"),
                     sep="\t", index=False)


@transform((buildWindowsFoldChangesPerMedian, buildWindowsFoldChangesPerInput),
           suffix(".tsv.gz"), ".load")
def loadWindowsFoldChanges(infile, outfile):
    '''
    Load fold change stats compared to input and median to database tables -
    <track>_l2foldchange_input and <track>_l2foldchange_median, where
    track is the prefix of the input bam filename.

    Parameters
    ----------
    infile: str
        name of :term:`tsv` file with fold change information
    outfile: str
        filename of database load log file.
    '''
    P.load(infile, outfile)


@transform((aggregateWindowsTagCounts, aggregateContextTagCounts),
           suffix(".tsv.gz"),
           "_stats.tsv")
def summarizeAllWindowsTagCounts(infile, outfile):
    '''
    Perform summarization of all read counts

    Parameters
    ----------
    infile: str
        filenames of :term:`tsv` formatted window and context based tag counts
    outfile: str
        filename of :term:`tsv` formatted file to write summarised tag counts
    '''

    prefix = P.snip(outfile, ".tsv")
    job_memory = "32G"

    statement = '''cgat runExpression
    --method=summary
    --tags-tsv-file=%(infile)s
    --output-filename-pattern=%(prefix)s_
    --log=%(outfile)s.log
    > %(outfile)s'''
    P.run()


@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"counts.dir/\1_stats.tsv")
def summarizeWindowsTagCounts(infiles, outfile):
    '''
    Perform summarization of read counts within experiments.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        design file in :term:`tsv` format
    infiles[1]: str
        aggregated tag counts for windows in :term:`tsv` format
    outfile: str
        filename of output file for summary of read counts in each experiment
    '''

    design_file, counts_file = infiles
    prefix = P.snip(outfile, ".tsv")
    statement = '''cgat runExpression
    --method=summary
    --design-tsv-file=%(design_file)s
    --tags-tsv-file=%(counts_file)s
    --output-filename-pattern=%(prefix)s_
    --log=%(outfile)s.log
    > %(outfile)s'''
    P.run()


@follows(mkdir("dump.dir"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"dump.dir/\1.tsv.gz")
def dumpWindowsTagCounts(infiles, outfile):
    '''
    Output tag tables used for debugging purposes. The tables
    can be loaded into R for manual analysis.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        design file in :term:`tsv` format
    infiles[1]: str
        aggregated tag counts for windows in :term:`tsv` format
    outfile: str
        :term:`tsv` filename for tag table dump

    '''
    design_file, counts_file = infiles

    statement = '''cgat runExpression
              --method=dump
              --design-tsv-file=%(design_file)s
              --tags-tsv-file=%(counts_file)s
              --log=%(outfile)s.log
              > %(outfile)s'''

    P.run()


@transform((summarizeWindowsTagCounts, summarizeAllWindowsTagCounts),
           suffix("_stats.tsv"), "_stats.load")
def loadTagCountSummary(infile, outfile):
    '''
    Load summaries of tag counts across windows across experiments (tables
    context_stats and windows_stats)
    and within experiments (context_count_stats and windows_count_stats) into
    the database.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`tsv` files contianing summarised tag
        counts
    outfile: str
        filename of database load log file

    '''
    P.load(infile, outfile)
    P.load(P.snip(infile, ".tsv") + "_correlation.tsv",
           P.snip(outfile, "_stats.load") + "_correlation.load",
           options="--first-column=track")


# @P.add_doc(PipelineWindows.normalizeBed)
@follows(buildWindows, countTagsWithinWindows)
@transform((aggregateWindowsTagCounts,
            aggregateContextTagCounts),
           suffix(".tsv.gz"),
           ".norm.tsv.gz")
def normalizeBed(infile, outfile):
    '''
    Normalize counts in a bed file by total library size.
    Return as bedGraph format

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` file containing tag counts within windows
    outfile: str
        filename of :term:`bedGraph` file to write results
    '''

    # normalize count column by total library size

    tmpfile = P.getTempFilename(shared=True)

    P.submit(module='CGATPipelines.PipelineWindows',
             function='normalizeBed',
             infiles=infile,
             outfiles=tmpfile,
             to_cluster=True,
             job_options="-l mem_free=32G")

    statement = '''cat %(tmpfile)s |
                   gzip > %(outfile)s; rm -f %(tmpfile)s'''

    P.run()


# @P.add_doc(PipelineWindows.enrichmentVsInput)
@follows(normalizeBed)
@transform("counts.dir/*.norm.bedGraph.gz",
           regex("counts.dir/(.+)-(.+)-(.+)_Input.bwa.norm.bedGraph.gz"),
           add_inputs(r"counts.dir/\1-\2-\3.bwa.norm.bedGraph.gz"),
           r"counts.dir/\1-\2-\3.vsInput.bedGraph.gz")
def enrichVsInput(infile, outfile):
    '''
    Calculate enrichment vs Input and output as bedGraph format

    Parameters
    ----------
    infile: list
        list of filenames
    infile[0]: str
        filename of normalised :term:`bedGraph` file showing counts in
        the input
    infile[1]: str
        filename of normalised :term:`bedGraph` files showing
        counts in each experiment
    outfile: str
        filename of output :term:`bedGraph` file
    '''

    tmpfile = P.getTempFilename(shared=True)
    P.submit(module='CGATPipelines.PipelineWindows',
             function='enrichmentVsInput',
             infiles=infile,
             outfiles=tmpfile,
             to_cluster=True)

    statement = '''cat %(tmpfile)s |  gzip > %(outfile)s; rm -f %(tmpfile)s'''

    P.run()


@follows(mkdir("bigwig.dir"), normalizeBed)
@transform("counts.dir/*.bedGraph.gz",
           regex("counts.dir/(.+).bedGraph.gz"),
           r"bigwig.dir/\1.bw")
def convertBed2BigWig(infile, outfile):
    '''
    Use UCSC tools to convert bedGraph -> bigwig

    Parameters
    ----------
    infile: str
        filename of :term:`bedGraph` formatted file
    outfile: str
        filename of :term:`bigwig` formatted file
    '''

    tmpfile = P.getTempFilename()

    contig_file = PARAMS['annotations_dir'] + "/contigs.tsv"

    statement = '''zcat %(infile)s | sort -k 1,1 -k 2,2n > %(tmpfile)s;
                   bedGraphToBigWig %(tmpfile)s %(contig_file)s %(outfile)s;
                   checkpoint ;
                   rm -f %(tmpfile)s'''

    P.run()


@follows(mkdir("images.dir"), convertBed2BigWig)
@transform(convertBed2BigWig,
           regex("bigwig.dir/(.+)-(.+)-(.+).bw"),
           r"images.dir/\1-\2-\3.hilbert.sentinel")
def plotHilbertCurves(infile, outfile):
    '''
    Use the BioC package `HilbertVis` to generate hilbert curves of bigwig
    files.  Generates one image file for each contig in the bigwig file.

    Parameters
    ----------
    infile: str
        filename of :term:`bigwig` formatted file
    outfile: str
        filename of hilbert plot
    '''
    statement = '''cgat bigwig2hilbert -v 0
                          --log=%(infile)s.log
                          --images-dir=images.dir
                          %(infile)s'''

    P.run()

    P.touch(outfile)


def loadMethylationData(infile, design_file):
    '''
    Load methylation data for deseq/edger analysis.

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups

    Parameters
    ----------
    infile: str
        filename of counts table in :term:`tsv` format to convert to an R
        object
    design_file: str
        filename of design file in :term:`tsv` format
    '''

    E.info("reading data")
    R('''counts_table=read.delim('%(infile)s', header=TRUE,'''
      '''row.names=1, stringsAsFactors=TRUE )''' % locals())

    E.info("read data: %i observations for %i samples" %
           tuple(R('''dim(counts_table)''')))

    # Load comparisons from file
    R('''pheno = read.delim('%(design_file)s', '''
      '''header = TRUE, stringsAsFactors = TRUE )''' % locals())

    # Make sample names R-like - substitute - for . and add the .prep suffix
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')

    # Ensure pheno rows match count columns
    R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),drop=FALSE]''')

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''conds <- pheno2$group[ includedSamples ]''')

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- factor(pheno2$group[ includedSamples ])''')
    R('''pairs = factor(pheno2$pair[ includedSamples ])''')

    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')

    E.info("filtered data: %i observations for %i samples" %
           tuple(R('''dim(countsTable)''')))

    return groups, pairs


# @P.add_doc(PipelineWindows.runDE)
@follows(mkdir("deseq.dir"), mkdir("deseq.dir/plots"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"deseq.dir/\1.tsv.gz")
def runDESeq(infiles, outfile):
    '''
    Estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.

    Parameters
    ----------
    infiles: list
    infiles[0]: str
        filename of design file in :term:`tsv` format
    infiles[1]: str
        filename of window tag count data in :term:`tsv` format

    deseq_fit_type: str
        :term:`PARAMS`
        fit type to estimate dispersion with deseq, refer to
        https://bioconductor.org/packages/release/bioc/manuals/DESeq/man/DESeq.pdf
    deseq_dispersion_method: str
        :term:`PARAMS`
        method to estimate dispersion with deseq, refer to
        https://bioconductor.org/packages/release/bioc/manuals/DESeq/man/DESeq.pdf
    deseq_sharing_mode: str
        :term:`PARAMS`
        determines which dispersion value is saved for each gene, refer to
        https://bioconductor.org/packages/release/bioc/manuals/DESeq/man/DESeq.pdf
    tags_filter_min_counts_per_row: int
        :term:`PARAMS`
        minimum number of counts below which to filter rows
    tags_filter_min_counts_per_sample: int
        :term:`PARAMS`
        minimum number of counts below which to filter samples
    tags_filter_percentile_rowsums: int
        :term:`PARAMS`
        percentile filtering using the total number of counts per row, e.g.
        20 removes 20% of windows with lowest counts.
    outfile: str
        filename of table to write deseq results in :term:`tsv` format
    '''

    spike_file = os.path.join("spike.dir", infiles[0]) + ".gz"
#    deseq_version = PARAMS['deseq_version']
    if os.path.exists(spike_file):
        outfile_spike = P.snip(outfile, '.tsv.gz') + '.spike.gz'
        PipelineWindows.runDE(infiles[0],
                              infiles[1],
                              outfile_spike,
                              "deseq.dir",
                              method="deseq",
                              spike_file=spike_file)
    PipelineWindows.runDE(infiles[0],
                          infiles[1],
                          outfile,
                          "deseq.dir",
                          method="deseq")


@transform(runDESeq, suffix(".tsv.gz"), ".load")
def loadDESeq(infile, outfile):
    '''Load DESeq per-chunk summary stats into database table <track>.

    Parameters
    ----------
    infile: str
        filename of Deseq output in :term:`tsv` format
    outfile: str
        logfile of database load
    '''

    prefix = P.snip(outfile, ".load")

    if os.path.exists(infile + "_size_factors.tsv"):
        P.load(infile + "_size_factors.tsv",
               prefix + "_deseq_size_factors.load",
               collapse=True,
               transpose="sample")

    for fn in glob.glob(infile + "*_summary.tsv"):
        prefix = P.snip(fn[len(infile) + 1:], "_summary.tsv")

        P.load(fn,
               prefix + ".deseq_summary.load",
               collapse=0,
               transpose="sample")

    P.touch(outfile)


# @P.add_doc(PipelineWindows.runDE)
@follows(mkdir("deseq2.dir"), mkdir("deseq2.dir/plots"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"deseq2.dir/\1.tsv.gz")
def runDESeq2(infiles, outfile):
    '''
    Estimate differential expression using DESeq2.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.

    Parameters
    ----------
    infiles: list
    infiles[0]: str
        filename of design file in :term:`tsv` format
    infiles[1]: str
        filename of window tag count data in :term:`tsv` format
    deseq2_fdr: float
        :term:`PARAMS`
        threshold fdr value for deseq2 analysis
    deseq2_model: str
        :term:`PARAMS`
        model to pass to deseq2, see
        https://bioconductor.org/packages/release/bioc/html/DESeq2.html
    deseq2_contrasts: str
        :term:`PARAMS`
        contrasts to return in pairwise tests in deseq2, see
        https://bioconductor.org/packages/release/bioc/html/DESeq2.html
    tags_filter_min_counts_per_row: int
        :term:`PARAMS`
        minimum number of counts below which to filter rows
    tags_filter_min_counts_per_sample: int
        :term:`PARAMS`
        minimum number of counts below which to filter samples
    tags_filter_percentile_rowsums: int
        :term:`PARAMS`
        percentile filtering using the total number of counts per row, e.g.
        20 removes 20% of windows with lowest counts.
    outfile: str
        filename of table to write deseq results in :term:`tsv` format
    '''

    spike_file = os.path.join("spike.dir", infiles[0]) + ".gz"

    if os.path.exists(spike_file):
        outfile_spike = P.snip(outfile, '.tsv.gz') + '.spike.gz'
        PipelineWindows.runDE(infiles[0],
                              infiles[1],
                              outfile_spike,
                              "deseq2.dir",
                              method="deseq2",
                              spike_file=spike_file)
    PipelineWindows.runDE(infiles[0],
                          infiles[1],
                          outfile,
                          "deseq2.dir",
                          method="deseq2")


@transform(runDESeq2, suffix(".tsv.gz"), ".load")
def loadDESeq2(infile, outfile):
    '''Load DESeq per-chunk summary stats into database table <track>.

    Parameters
    ----------
    infile: str
        filename of Deseq2 output in :term:`tsv` format
    outfile: str
        logfile of database load
    '''

    prefix = P.snip(outfile, ".load")

    if os.path.exists(infile + "_size_factors.tsv"):
        P.load(infile + "_size_factors.tsv",
               prefix + "_deseq2_size_factors.load",
               collapse=True,
               transpose="sample")

    for fn in glob.glob(infile + "*_summary.tsv"):
        prefix = P.snip(fn[len(infile) + 1:], "_summary.tsv")

        P.load(fn,
               prefix + ".deseq2_summary.load",
               collapse=0,
               transpose="sample")

    P.touch(outfile)


@follows(mkdir("spike.dir"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"spike.dir/\1.tsv.gz")
def buildSpikeIns(infiles, outfile):
    '''Build a table with counts to spike into the original count
    data sets.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        design table in :term:`tsv` format
    infiles[1]: str
        filename of windows tag counts in :term:`tsv` format
    outfile: str
        filename of spike in file in :term:`tsv` format
    '''

    design_file, counts_file = infiles
    design = P.snip(design_file, ".tsv")
    statement = '''
    zcat %(counts_file)s
    | cgat runExpression
            --log=%(outfile)s.log
            --design-tsv-file=%(design_file)s
            --tags-tsv-file=-
            --method=spike
            --output-filename-pattern=%(outfile)s_
    | gzip
    > %(outfile)s
    '''
    P.run()


# @P.add_doc(PipelineWindows.runDE)
@follows(mkdir("edger.dir"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"edger.dir/\1.tsv.gz")
def runEdgeR(infiles, outfile):
    '''estimate differential methylation using EdgeR

    This method applies a paired test. The analysis follows
    the example in chapter 11 of the EdgeR manual.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        design table in :term:`tsv` format
    infiles[1]: str
        filename of windows tag counts in :term:`tsv` format
    edger_dispersion: float
        :term:`PARAMS`
        typical dispersion when there are no replicates, refer to
        https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    edger_fdr: float
        :term:`PARAMS`
        false discovery rate threshold, refer to
        https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    tags_filter_min_counts_per_row: int
        :term:`PARAMS`
        minimum number of counts below which to filter rows
    tags_filter_min_counts_per_sample: int
        :term:`PARAMS`
        minimum number of counts below which to filter samples
    tags_filter_percentile_rowsums: int
        :term:`PARAMS`
        percentile filtering using the total number of counts per row, e.g.
        20 removes 20% of windows with lowest counts.
    outfile: str
        filename of edger output file in :term:`tsv` format
    '''

    spike_file = os.path.join("spike.dir", infiles[0]) + ".gz"
    if os.path.exists(spike_file):
        outfile_spike = P.snip(outfile, '.tsv.gz') + '.spike.gz'

        PipelineWindows.runDE(infiles[0],
                              infiles[1],
                              outfile_spike,
                              "edger.dir",
                              method="edger",
                              spike_file=spike_file)

    PipelineWindows.runDE(infiles[0],
                          infiles[1],
                          outfile,
                          "edger.dir",
                          method="edger")


@transform(runEdgeR, suffix(".tsv.gz"), ".load")
def loadEdgeR(infile, outfile):
    '''load EdgeR per-chunk summary stats

    Parameters
    ----------
    infile: str
        filename of edger output in :term:`tsv` format
    outfile: str
        logfile of database load
    '''

    prefix = P.snip(outfile, ".load")

    for fn in glob.glob(infile + "*_summary.tsv"):
        prefix = P.snip(fn[len(infile) + 1:], "_summary.tsv")

        P.load(fn,
               prefix + ".edger_summary.load",
               collapse=0,
               transpose="sample")

    P.touch(outfile)


# @P.add_doc(PipelineWindows.outputRegionsOfInterest)
@follows(mkdir("roi.dir"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           add_inputs(aggregateWindowsTagCounts),
           r"roi.dir/\1.tsv.gz")
def runFilterAnalysis(infiles, outfile):
    '''
    Output windows applying a filtering criterion.
    Does not apply a threshold . (?)

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        design table in :term:`tsv` format
    infiles[1]: str
        filename of windows tag counts in :term:`tsv` format
    outfile: str
        filename of filtering output file in :term:`tsv` format
    '''
    PipelineWindows.outputRegionsOfInterest(
        infiles[0],
        infiles[1],
        outfile)


# @P.add_doc(PipelineWindows.runMEDIPSDMR)
@follows(mkdir("medips.dir"))
@transform("design*.tsv",
           regex("(.*).tsv"),
           r"medips.dir/\1.tsv.gz")
def runMedipsDMR(infile, outfile):
    '''
    Run MEDIPS single file analysis

    Parameters
    ----------
    infile: str
        filename of design file in :term:`tsv` format
    medips_genome: str
        :term:`PARAMS`
        UCSC genome name using R naming convention e.g. Hsapiens.UCSC.hg19
    medips_shift: str
        :term:`PARAMS`
        shift between windows to use in medips.
        Tags represent the ends of fragments - shift towards the
        3' direction to improve their representation of the binding site.
        Tools are available to calculate this.
    medips_extension: str
        :term:`PARAMS`
         extend reads to the length of sonicated fragments
    medips_window_size: int
        :term:`PARAMS`
        window size to use for MEDIPS analysis
    medips_fdr: float
        :term:`PARAMS`
        threshold for false discovery rate for MEDIPS
    outfile: str
        filename of :term:`tsv` formatted file to write the medips output
        table


    '''
    PipelineWindows.runMEDIPSDMR(infile, outfile)


DIFFTARGETS = []
mapToTargets = {'deseq': (loadDESeq, runDESeq,),
                'edger': (runEdgeR,),
                'filter': (runFilterAnalysis,),
                'medips': (runMedipsDMR,),
                'deseq2': (loadDESeq2, runDESeq2)
                }
for x in METHODS:
    DIFFTARGETS.extend(mapToTargets[x])


@follows(loadTagCountSummary,
         loadWindowStats,
         loadWindowsTagCounts,
         loadWindowsFoldChanges,
         *DIFFTARGETS)
def diff_windows():
    '''
    Records when all differential expression analysis is complete
    '''
    pass


@transform(DIFFTARGETS, suffix(".gz"), ".cpg.tsv.gz")
def computeWindowComposition(infile, outfile):
    '''
    For the windows returned from differential analysis, compute CpG
    content for QC purposes.

    Parameters
    ----------
    infile: str
        output from differential expression with DEseq, EdgeR, filtering
        and/or MEDIPS.
    outfile: str
        filename of :term:`tsv` formatted file to write cpg composition data.
    '''

    statement = '''
    zcat %(infile)s
    | grep -v "^spike"
    | perl -p -e "s/:/\\t/; s/-/\\t/; s/test_id/contig\\tstart\\tend/"
    | cgat bed2table
    --counter=composition-cpg
    --genome-file=%(genome_dir)s/%(genome)s
    --has-header
    | gzip
    > %(outfile)s
    '''

    P.run()


@transform(computeWindowComposition, suffix(".tsv.gz"), ".load")
def loadWindowComposition(infile, outfile):
    '''
    Load a sample of window composition data for QC purposes in database
    table <track_cpg> where track is the name of the differential expression
    output file

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` formatted file containing cpg composition data.
    outfile: str
        filename of database load logfile
    '''
    P.load(infile, outfile, limit=10000)


@transform(DIFFTARGETS,
           suffix(".tsv.gz"),
           ".linear")
def outputGWASFiles(infile, outfile):
    '''
    Output GWAS formatted files for viewing in the
    IGV genome browser. The output files contain
    chromosomal position and an associated P-Value.

    See here for acceptable file formats:
    http://www.broadinstitute.org/software/igv/GWAS

    Parameters
    ----------
    infile: str
        output from differential expression with DEseq, EdgeR, filtering
        and/or MEDIPS.
    outfile: str
        filename of :term:`GWAS` formatted file to write data.
    '''

    statement = '''
    zcat %(infile)s
    | grep -v "^spike"
    | awk 'BEGIN {printf("chr\\tpos\\tsnp\\tpvalue\\n")}
    !/^test_id/ {split($1,a,":");
    printf("%%s\\t%%i\\t%%s\\t%%s\\n", a[1], a[2], $1, $8)}'
    | %(pipeline_scriptsdir)s/hsort 1 -k1,1 -k2,2n
    > %(outfile)s
    '''

    P.run()


@transform(DIFFTARGETS,
           suffix(".tsv.gz"),
           ".merged.tsv.gz")
def mergeDMRWindows(infile, outfile):
    '''
    Merge overlapping windows.

    Sample/control labels are by default inverted to reflect
    that unmethylated windows are of principal interest.

    Parameters
    ----------
    infile: str
        output from differential expression with DEseq, EdgeR, filtering
        and/or MEDIPS.
    outfile: str
        filename of :term:`tsv` formatted file to write merged window data.
    '''
    # the other outfiles will be created automatically by
    # the script medip_merge_intervals

    prefix = P.snip(outfile, ".tsv.gz")

    job_memory = "3G"

    statement = '''
    zcat %(infile)s
    | grep -v "^spike"
    | cgat medip_merge_intervals
    --log=%(outfile)s.log
    --invert
    --output-filename-pattern=%(prefix)s.%%s.bed.gz
    | gzip
    > %(outfile)s
    '''

    P.run()


# @P.add_doc(PipelineWindows.buildSpikeResults)
@transform(DIFFTARGETS, suffix(".gz"), ".power.gz")
def buildSpikeResults(infile, outfile):
    '''
    Compute results of spike analysis and put into a matrix.

    Parameters
    ----------
    infile: str
        output from differential expression with DEseq, EdgeR, filtering
        and/or MEDIPS.
    outfile: str
        output filename in :term:`tsv` format for spike in results
    '''
    PipelineWindows.buildSpikeResults(infile, outfile)


@transform(buildSpikeResults, suffix('.tsv.power.gz'), '.power.load')
def loadSpikeResults(infile, outfile):
    '''load spike in results to database table <track>.power where track
    is differential expression output prefix.

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` formatted file containing spike in data
    outfile: str
        filename of database load logfile
    '''
    method = P.snip(os.path.dirname(outfile), '.dir')
    tablename = P.toTable(outfile)
    tablename = '_'.join((tablename, method))

    P.load(infile, outfile, options='--add-index=fdr,power --allow-empty-file',
           tablename=tablename)


# @P.add_doc(PipelineWindows.buildDMRStats)
@transform(mergeDMRWindows, suffix(".merged.tsv.gz"), ".stats")
def buildDMRStats(infile, outfile):
    '''Compute differential methylation stats - count the number of up/down,
    2fold up/down etc. results

    Parameters
    ----------
    infile: str
        filename in :term:`tsv` format  with differential expression output
    outfile: str
        :term:`tsv`filename to write differential methylation stats

    '''
    method = os.path.dirname(infile)
    method = P.snip(method, ".dir")
    PipelineWindows.buildDMRStats([infile], outfile, method=method)


# @P.add_doc(PipelineWindows.plotDETagStats)
@transform(mergeDMRWindows,
           suffix(".tsv.gz"),
           add_inputs(buildWindowComposition),
           ".plots")
def plotDETagStats(infiles, outfile):
    '''Plot differential expression stats

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        filename of DMR analysis output in :term:`tsv` format
    infiles[1]: str
        filename of window composition table in :term:`tsv` format
    outfile: str
        output filename, dummy filename used to trigger next step in pipeline.
    '''

    PipelineWindows.plotDETagStats(
        infiles[0], infiles[1], outfile,
        submit=True,
        job_options="-l mem_free=16")


# @P.add_doc(PipelineWindows.buildFDRStats)
@transform(mergeDMRWindows, suffix(".merged.tsv.gz"), ".fdr")
def buildFDRStats(infile, outfile):
    '''Compute differential methylation false discovery rate stats

    Parameters
    ----------
    infile: str
        filename of DMR analysis output in :term:`tsv` format
    outfile: str
        :term:`tsv`filename to write FDR stats

    '''
    method = os.path.dirname(infile)
    method = P.snip(method, ".dir")
    PipelineWindows.buildFDRStats(infile, outfile, method=method)


# @P.add_doc(PipelineWindows.outputAllWindows)
@transform(mergeDMRWindows,
           suffix(".merged.tsv.gz"),
           ".mergedwindows.all.bed.gz")
def outputAllWindows(infile, outfile):
    '''output all windows as a bed file with the l2fold change as a score

    Parameters
    ----------
    infile : string
       Input filename in :term:`tsv` format. Typically the output
       from :mod:`scripts/runExpression`.
    outfile : string
       Output filename in :term:`bed` format.

    '''
    PipelineWindows.outputAllWindows(infile, outfile)


@transform(outputAllWindows, suffix(".all.bed.gz"),
           (".top.bed.gz", ".bottom.bed.gz"))
def outputTopWindows(infile, outfiles):
    '''
    Output :term:`bed` file with largest/smallest l2fold changes.
    The resultant bed files are sorted by coordinate.

    Parameters
    ----------
    infile: str
        filename of :term:`bed` file containing windows and l2fold changes
    outfiles: list
        list of filenames
    outfiles[0]: str
        filename of :term:`bed` file to write largest l2fold changes
    outfiles[1]: str
        filename of :term:`bed` file to write smallest l2fold changes
    '''
    outfile = outfiles[0]

    ignore_pipe_errors = True

    statement = '''zcat %(infile)s
    | awk '$4 !~ /inf/'
    | sort -k4,4n
    | tail -n %(bed_export)i
    | sort -k1,1 -k2,2n
    | bgzip
    > %(outfile)s
    '''
    P.run()

    outfile = outfiles[1]

    statement = '''zcat %(infile)s
    | awk '$4 !~ /inf/'
    | sort -k4,4n
    | head -n %(bed_export)i
    | sort -k1,1 -k2,2n
    | bgzip
    > %(outfile)s
    '''
    P.run()


@transform(mergeDMRWindows,
           suffix(".merged.tsv.gz"),
           ".stats")
def buildDMRWindowStats(infile, outfile):
    '''Compute window size statistics of DMR from bed file.

    Parameters
    ----------
    infile: str
        filename of DMR analysis output in :term:`tsv` format
    outfile: str
        filename in :term:`tsv` format to write window size statistics as
        histograms
    '''

    statement = '''
    zcat %(infile)s
    | grep -v 'contig'
    | cgat gff2histogram
                   --force-output
                   --format=bed
                   --output-section=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()


# @P.add_doc(PipelineWindows.buildDMRStats)
@transform(runMedipsDMR, suffix(".tsv.gz"), ".stats")
def buildMedipsStats(infile, outfile):
    '''Compute differential methylation stats for medips data

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` format file containing MEDIPs output
    medips_fdr: float
        :term:`PARAMS`
        threshold false discovery rate for medips
    outfile: str
        filename of :term:`tsv` file to write medips DMR stats
    '''
    method = os.path.dirname(infile)
    method = P.snip(method, ".dir")
    infiles = glob.glob(infile + "*_data.tsv.gz")
    PipelineWindows.buildDMRStats(infiles,
                                  outfile,
                                  method=method,
                                  fdr_threshold=PARAMS["medips_fdr"])


@merge(buildDMRStats, "dmr_stats.load")
def loadDMRStats(infiles, outfile):
    '''load DMR stats into database table - dmr_stats

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` formatted file containing DMR stats
    outfile: str
        filename of database load logfile
    '''
    P.concatenateAndLoad(infiles, outfile,
                         missing_value=0,
                         regex_filename=".*\/(.*).stats")

# @merge( buildDMRBed, "dmr_overlap.tsv.gz" )
# def computeDMROverlap( infiles, outfile ):
#     '''compute overlap between bed sets.'''

#     to_cluster = True

#     if os.path.exists(outfile):
# note: update does not work due to quoting
#         os.rename( outfile, "orig." + outfile )
#         options = "--update=orig.%s" % outfile
#     else:
#         options = ""

#     infiles = " ".join( infiles )

# note: need to quote track names
#     statement = '''
#         cgat diff_bed
#               --pattern-identifier=".*/(.*).dmr.bed.gz"
#               --log=%(outfile)s.log
#               %(options)s %(infiles)s
# | awk -v OFS="\\t" '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
#         | gzip
#         > %(outfile)s
#         '''

#     P.run()


@transform(mergeDMRWindows, regex("(.*)\.(.*).merged.gz"), r"\1_\2.bed.gz")
def buildMRBed(infile, outfile):
    '''output bed6 file with methylated regions.
    All regions are output, even the insignificant ones.
    The score is the log fold change.

    Parameters
    ----------
    infile: str
        filename of DMR analysis output in :term:`tsv` format
    outfile: str
        filename of :term:`bed6` file to write methylated regions
    '''

    outf = IOTools.openFile(outfile, "w")
    c = E.Counter()
    for row in csv.DictReader(IOTools.openFile(infile),
                              dialect="excel-tab"):
        c.input += 1

        contig, start, end = re.match(
            "(.*):(\d+)-(\d+)", row["interval_id"]).groups()
        c.output += 1
        outf.write(
            "\t".join((contig, start, end, str(c.input), row["lfold"])) + "\n")

    outf.close()

    E.info("%s" % str(c))


@follows(mkdir("overlaps.dir"), mergeDMRWindows)
@collate(("edger.dir/*merged*.bed.gz",
          "deseq.dir/*merged*.bed.gz",
          "medips.dir/*merged*.bed.gz"),
         regex("(.*).dir/(.*).merged.(.*).bed.gz"),
         r"overlaps.dir/method_\2_\3.overlap")
def buildOverlapByMethod(infiles, outfile):
    '''Compute overlap between differential expression output using different
    methods.

    Parameters
    ----------
    infiles: list
        list of filenames of DMR output files in :term:`bed` format
    outfile: str
        filename to write overlap information

    '''

    if os.path.exists(outfile):
        # note: update does not work due to quoting
        os.rename(outfile, outfile + ".orig")
        options = "--update=%s.orig" % outfile
    else:
        options = ""

    infiles = " ".join(infiles)

    # note: need to quote track names
    statement = '''
    cgat diff_bed %(options)s %(infiles)s
    | awk -v OFS="\\t"
    '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
    > %(outfile)s
    '''

    P.run()


@follows(mkdir("overlaps.dir"), mergeDMRWindows)
@collate(("edger.dir/*merged*.bed.gz",
          "deseq.dir/*merged*.bed.gz",
          "medips.dir/*merged*.bed.gz"),
         regex("(.*).dir/(.*).merged.(.*).bed.gz"),
         r"overlaps.dir/\1_\3.overlap")
def buildOverlapWithinMethod(infiles, outfile):
    '''Compute overlap between differential expression output for different
    files using the same method.

    Parameters
    ----------
    infiles: list
        list of filenames of DMR output files in :term:`bed` format
    outfile: str
        filename to write overlap information

    '''
    if os.path.exists(outfile):
        # note: update does not work due to quoting
        os.rename(outfile, outfile + ".orig")
        options = "--update=%s.orig" % outfile
    else:
        options = ""

    infiles = " ".join(infiles)

    # note: need to quote track names
    statement = '''
    cgat diff_bed %(options)s %(infiles)s
    | awk -v OFS="\\t"
    '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
    > %(outfile)s
    '''

    P.run()


@transform((buildOverlapByMethod,
            buildOverlapWithinMethod),
           suffix(".overlap"), "_overlap.load")
def loadOverlap(infile, outfile):
    '''Load results of overlap analyses within and between methods to database
    table <track>_overlap where track is the prefix of the overlap analysis
    filename

    Parameters
    ----------
    infile: str
        filename containing the results of an overlap analysis
    outfile: str
        database load logfile
    '''
    P.load(infile, outfile,
           tablename="overlap",
           options="--add-index=set1 --add-index=set2")


@follows(loadDMRStats, loadSpikeResults,
         outputAllWindows, outputTopWindows)
def dmr():
    '''
    Records when all DMR analysis in pipeline is complete.
    '''
    pass


@follows(mergeDMRWindows)
@merge(("edger.dir/*.merged.*.bed.gz",
        "deseq.dir/*.merged.*.bed.gz",
        "medips.dir/*.merged.*.bed.gz"),
       "dmr.bed.gz")
def combineWindows(infiles, outfile):
    '''
    Build a :term:`bed` file that contains all intervals
    detected as DMR windows. The name field will
    be set to the filename.

    Parameters
    ----------
    infiles: list
        list of filenames of DMR analysis results in :term:`bed` format
    outfile: str
        filename of :term:`bed` file to contain all DMR windows
    '''

    statement = []
    if os.path.exists(outfile):
        os.unlink(outfile)
    for x in infiles:
        statement.append("""zcat %(x)s
        | awk '{printf("%%s\\t%%i\\t%%i\\t%(x)s")}'
        | gzip
        >> %(outfile)s""" % locals())
    statement = ";".join(statement)
    P.run(statement)


# @P.add_doc(PipelineWindows.summarizeTagsWithinContext)
@follows(mkdir('contextstats.dir'))
@transform(prepareTags,
           regex(".*/(.*).bed.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_genomic_context_bed"])),
           r"contextstats.dir/\1.contextstats.tsv.gz")
def summarizeTagsWithinContext(infiles, outfile):
    """
    Count number of tags overlapping various genomic contexts (by 50% of more).
    Counts are summarized for each category.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        filename of :term:`bed` file containing tag counts
    infiles[1]: str
        filename of :term:`bed` file containing genomic context
    outfile: str
        filename of :term:`tsv` file to write the context stats
    """
    PipelineWindows.summarizeTagsWithinContext(
        infiles[0],
        infiles[1],
        outfile)


# @P.add_doc(PipelineWindows.mergeSummarizedContextStats)
@merge(summarizeTagsWithinContext, "contextstats.dir/contextstats.tsv.gz")
def mergeSummarizedContextStats(infiles, outfile):
    """Merge information on tags overlapping various genomic contexts from
    different samples into :term:`tsv` file.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`tsv` formatted information about tags
        overlapping genomic contexts for each sample
    outfile: str
        filename for :term:`tsv` formatted file to write merged context data
    """
    PipelineWindows.mergeSummarizedContextStats(infiles, outfile)


# @P.add_doc(PipelineWindows.loadSummarizedContextStats)
@merge(summarizeTagsWithinContext, "context_stats.load")
def loadSummarizedContextStats(infiles, outfile):
    """Load context mapping statistics into database table context_stats

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`tsv` formatted information about tags
        overlapping genomic contexts for each sample
    outfile: str
        database load logfile
    """
    PipelineWindows.loadSummarizedContextStats(infiles, outfile)


# @P.add_doc(PipelineWindows.mergeSummarizedContextStats)
@merge(summarizeTagsWithinContext,
       "contextstats.dir/contextstats_counts.tsv.gz")
def mergeSummarizedContextStatsAsCounts(infiles, outfile):
    '''Merge information on tags overlapping various genomic contexts from
    different samples into :term:`tsv` file with samples as columns rather
    than rows ready for normalisation.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`tsv` formatted information about tags
        overlapping genomic contexts for each sample
    outfile: str
        filename for :term:`tsv` formatted file to write merged context data
    '''
    PipelineWindows.mergeSummarizedContextStats(infiles, outfile,
                                                samples_in_columns=True)


# @P.add_doc(PipelineWindows.normalizeTagCounts)
@transform(mergeSummarizedContextStatsAsCounts,
           suffix(".tsv.gz"),
           "_normed.tsv.gz")
def normalizeSummarizedContextStats(infile, outfile):
    """Normalize context mapping statistics by dividing by the row-wise total
    (total across all samples)

    Parameters
    ----------
    infile: str
        filename of :term:`tsv` file with one sample per column summarising
        information about tags overlapping different genomic contexts
    outfile: str
        filename of :term:`tsv` file to write the normalised context
        statistics.
    """

    PipelineWindows.normalizeTagCounts(
        infile,
        outfile,
        method="total-row")


@follows(mkdir("gat_context.dir"))
@transform(prepareTags,
           regex(".*/(.*).bed.gz"),
           add_inputs(
               os.path.join(
                   PARAMS["annotations_dir"],
                   PARAMS["annotations_interface_genomic_context_bed"]),
               os.path.join(
                   PARAMS["annotations_dir"],
                   PARAMS["annotations_interface_contigs_ungapped_bed"])),
           r"gat_context.dir/\1.tsv.gz")
def testTagContextOverlap(infiles, outfile):
    """test for genomic overlap using gat."""

    tagfile, contextfile, workspacefile = infiles
    PipelineWindows.testTagContextOverlap(
        tagfile,
        contextfile,
        workspacefile,
        outfile,
        job_threads=PARAMS["gat_threads"],
        samples=PARAMS["gat_samples"],
        options=PARAMS["gat_options"])


@transform(testTagContextOverlap,
           suffix(".tsv.gz"),
           "_context_gat.load")
def loadTagContextOverlap(infile, outfile):
    """load results from GAT analysis into database.
    """
    P.load(infile, outfile, options="--add-index=track")


@follows(mkdir("medips.dir"))
@transform("*.bam",
           regex("(.*).bam"),
           r"medips.dir/\1.tsv.gz")
def runMedipsQC(infile, outfile):
    '''run MEDIPS single file analysis
    '''
    PipelineWindows.runMEDIPSQC(infile, outfile)


@follows(mkdir("transcriptprofiles.dir"))
@transform(prepareTags,
           regex(r".*/([^/].*)\.bed.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_coding_exons_gtf"])),
           r"transcriptprofiles.dir/\1.transcriptprofile.tsv.gz")
def buildTranscriptProfiles(infiles, outfile):
    '''
    Build a table with the overlap profile of tags with protein coding exons.

    Parameters
    ----------
    infiles: list
        list of filenames
    infiles[0]: str
        filename of tag counts in :term:`tsv` format
    infiles[1]: str
        filename of annotation of exons in :term:`gtf` format
    outfile: str
        filename to write the overlap profile in :term:`tsv` format
    '''

    bedfile, gtffile = infiles

    track = P.snip(os.path.basename(outfile), '.transcriptprofile.tsv.gz')
    try:
        t = Sample(filename=track)
    except ValueError as msg:
        print(msg)
        return

    job_memory = "8G"

    # no input normalization, this is done later.
    options = ''
    # input_files = getInput( t )

    # currently only implement one input file per track
    # assert len(input_files) <= 1,\
    # "%s more than input: %s" % (track, input_files)

    # if len(input_files) == 1:
    #     options = '--controlfile=%s' % \
    #         (os.path.join( os.path.dirname( bedfile ),
    #                        input_files[0] + '.bed.gz') )
    statement = '''zcat %(gtffile)s
                   | cgat gtf2gtf
                     --method=filter
                     --filter-method=representative-transcript
                     --log=%(outfile)s.log
                   | cgat bam2geneprofile
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --reporter=transcript
                      --method=geneprofile
                      --method=tssprofile
                      --method=separateexonprofilewithintrons
                      --normalize-profile=all
                      --output-all-profiles
                      --resolution-upstream=1000
                      --resolution-downstream=1000
                      --resolution-cds=1000
                      --resolution-first-exon=1000
                      --resolution-last-exon=1000
                      --resolution-introns=1000
                      --extension-upstream=5000
                      --extension-downstream=5000
                      %(options)s
                      %(bedfile)s -
                   > %(outfile)s
                '''
    P.run()


@follows(loadTagContextOverlap, loadSummarizedContextStats)
def context():
    """Indicates that the context_stats part of the pipeline is complete"""
    pass


@follows(buildTranscriptProfiles,
         gc,
         loadPicardDuplicateStats,
         diff_windows,
         dmr)
def full():
    pass


@follows(buildBackgroundWindows)
def test():
    pass


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
    '''publish files.'''

    # directory : files
    export_files = {
        "bedfiles":
        glob.glob("deseq.dir/*.bed.gz") +
        glob.glob("edger.dir/*.bed.gz"),
    }

    # publish web pages
    P.publish_report(export_files=export_files)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
