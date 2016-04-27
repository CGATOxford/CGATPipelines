"""
pipeline_peakcalling.py - Window based genomic analysis
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python


Methods
=======

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----


Pipeline output
===============


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
import numpy
import sqlite3
import pandas
import glob
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelinePeakcallingScrum as PipelinePeakcalling
import CGAT.BamTools as Bamtools

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

BAMS = ['*.bam']

for fname in os.listdir("."):
    if fname.endswith('bam'):
        testbam = fname
        break
if Bamtools.isPaired(testbam) is True:
    PARAMS['paired_end'] = True

###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks

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

###############################################################
# Preprocessing Steps
@follows(mkdir("filtered_bams.dir"))
@transform(BAMS, regex("(.*).bam"), [r"filtered_bams.dir/\1_filtered.bam",
                                     r"filtered_bams.dir/\1_counts.tsv"])
def filterBAM(infile, outfiles):
    '''
    Applies various filters specified in the pipeline.ini to the bam file
    Currently implemented are filtering unmapped, unpaired and duplicate
    reads and filtering reads overlapping with blacklisted regions
    based on a list of bam files.
    '''
    filters = PARAMS['filters_bamfilters'].split(",")
    bedfiles = PARAMS['filters_bedfiles'].split(",")
    blthresh = PARAMS['filters_blacklistthresh']
    PipelinePeakcalling.filterBams(infile, outfiles, filters, bedfiles,
                                   float(blthresh),
                                   PARAMS['paired_end'],
                                   PARAMS['filters_strip'],
                                   PARAMS['filters_keepint'])


@transform(filterBAM, suffix("_filtered.bam"), "_insertsize.tsv")
def estimateInsertSize(infiles, outfile):
    '''
    Predicts insert size using MACS2 for single end data and using Bamtools
    for paired end data.
    Output is stored in insert_size.tsv
    '''
    infile = infiles[0]
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


if int(PARAMS['IDR_run']) == 1:
    @follows(mkdir("pseudo_bams.dir"))
    @subdivide(filterBAM, regex("filtered_bams.dir/(.*).bam"),
               [r"pseudo_bams.dir/\1_pseudo_1.bam",
                r"pseudo_bams.dir/\1_pseudo_2.bam"])
    def makePseudoBams(infiles, outfiles):
        '''
        Generates pseudo bam files each containing approximately 50% of reads
        from the original bam file for IDR analysis.

        '''
        infile = infiles[0]
        PipelinePeakcalling.makePseudoBams(infile, outfiles,
                                           PARAMS['paired_end'],
                                           PARAMS['IDR_randomseed'],
                                           submit=True)
else:
    @transform(filterBAM, regex("filtered_bams.dir/(.*).bam"),
               r'filtered_bams.dir/\1.bam')
    def makePseudoBams(infile, outfile):
        '''
        Dummy task if IDR not requested.
        '''
        pass


@follows(mergeInsertSizes)
@transform(makePseudoBams, regex("(.*)_bams\.dir\/(.*)\.bam"),
           r"\1_bams.dir/\2.bam")
def preprocessing(infile, outfile):
    '''
    Dummy task to ensure all preprocessing has run and
    bam files are passed individually to the next stage.
    '''
    pass


################################################################
# Peakcalling Steps

###############################################################
# The following are variable created just for testing purposes
###############################################################

test_inf = "./neural-SMC3-1.genome.bam"
control_inf = "./neural-input-1.genome.bam"
#control_inf = None
contigsfile = "/ifs/mirror/annotations/hg19_ensembl75_hierarchical/assembly.dir/contigs.tsv"


@transform(test_inf,
           regex("./neural-SMC3-1.genome.bam"),
           "./test.out")
def callMacs2peaks(infile, outfile):

    peakcaller = PipelinePeakcalling.Macs2Peakcaller(
        threads=1, paired_end=True, output_all=True,
        tool_options=PARAMS['macs2_options'], tagsize=None)
    statement = peakcaller.build(infile, outfile, contigsfile, control_inf)

    P.run()

################################################################
# QC Steps


################################################################

def full():
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

    # publish web pages
    P.publish_report(export_files=export_files)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
