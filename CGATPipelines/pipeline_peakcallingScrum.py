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
from ruffus import transform, merge, mkdir, follows, \
    regex, suffix, add_inputs, collate

from ruffus import *

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

import CGAT.BamTools as BamTools
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
@transform(BAMS, regex("(.*).bam"), r"filtered_bams.dir/\1")
def filter(infile, outfile):
    filters = PARAMS['filters_bamfilters']
    PipelinePeakcalling.filterBams(infile, outfile, filters)

################################################################
# Peakcalling Steps

###############################################################
# The following are variable created just for testing purposes
###############################################################

test_inf = "neural-SMC3-1.genome.bam"
control_inf = "neural-input-1.genome.bam"
# control_inf = None
contigsfile = ("/ifs/mirror/annotations/"
               "hg19_ensembl75_hierarchical/assembly.dir/contigs.tsv")


@transform(test_inf,
           regex("neural-SMC3-1.genome.bam"),
           "test.out.macs2")
def callMacs2peaks(infile, outfile):

    peakcaller = PipelinePeakcalling.Macs2Peakcaller(
        threads=1,
        paired_end=True,
        output_all=True,
        tool_options=PARAMS['macs2_options'],
        tagsize=None)
    statement = peakcaller.build(infile, outfile, contigsfile, control_inf)

    print statement
    P.run()


# TS - delete when Preprocess fragment size estimation had been performed
@follows(mkdir('fragment_size.dir'))
@transform(test_inf,
           regex("(\S+).genome.bam"),
           r"fragment_size.dir/\1.fragment_size")
def predictFragmentSize(infile, outfile):
    '''predict fragment size.

    For single end data, use MACS2, for paired-end data
    use the bamfile.

    In some BAM files the read length (rlen) attribute
    is not set, which causes MACS2 to predict a 0.

    Thus it is computed here from the CIGAR string if rlen is 0.
    '''
    tagsize = BamTools.estimateTagSize(infile, multiple="mean")

    if BamTools.isPaired(infile):
        mode = "PE"
        mean, std, n = BamTools.estimateInsertSizeDistribution(infile, 10000)
    else:
        mode = "SE"
        statement = '''macs2 predictd
        --format BAM
        --ifile %(infile)s
        --outdir %(outfile)s.dir
        --verbose 2
        %(macs2_predictd_options)s
        >& %(outfile)s.tsv
        '''
        P.run()

        with IOTools.openFile(outfile + ".tsv") as inf:
            lines = inf.readlines()
            line = [x for x in lines
                    if "# predicted fragment length is" in x]
            if len(line) == 0:
                raise ValueError(
                    'could not find predicted fragment length')
            mean = re.search(
                "# predicted fragment length is (\d+)",
                line[0]).groups()[0]
            std = 'na'

    outf = IOTools.openFile(outfile, "w")
    outf.write("mode\tfragmentsize_mean\tfragmentsize_std\ttagsize\n")
    outf.write("\t".join(
        map(str, (mode, mean, std, tagsize))) + "\n")
    outf.close()
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
