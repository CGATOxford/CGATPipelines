##########################################################################
#
#   MRC FGU Computational Genomics Analysis & Training Programme
#
#   $Id$
#
#   Copyright (C) 2014 David Sims
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
##########################################################################

"""====================
ReadQc pipeline
====================

:Author: David Sims
:Date: |today|
:Tags: Python

The readqc pipeline imports unmapped reads from one or more input
files and performs basic quality control steps. The pipeline performs
also read pre-processing.

Quality metrics are based on the fastqc tools, see
see http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ for further details.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

When pre-processing reads before mapping, the workflow of the pipeline is
as follows:

1. Run the ``full`` target to perform initial QC on the raw data. Then
   build the report (``build_report`` target).

2. Inspect the output to decide if and what kind of pre-processing is
   required.

3. Edit the configuration file ``pipeline.ini`` to activate
   pre-processing and parameterize it appropriately. Note that
   parameters can be set on a per-sample basis.

4. Rerun the ``full`` target and ``build_report`` targets. The data
   will now be processed and additional QC will be performed on the
   processed data. Note that all the processed data will be found in
   the :file:`processed.dir` directory.

Configuration
-------------

No general configuration required.

Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The default file format assumes the following convention:

   <sample>.<suffix>

The ``suffix`` determines the file type. The following suffixes/file
types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:
   `fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Pipeline output
----------------

The major output is a set of HTML pages and plots reporting on the quality of
the sequence archive

Example
=======

Example data is available at

https://www.cgat.org/downloads/public/cgatpipelines/pipeline_test_data/test_readqc.tgz

To run the example, simply unpack and untar::

   wget -qO- https://www.cgat.org/downloads/public/cgatpipelines/pipeline_test_data/test_readqc.tgz | tar -xvz
   cd test_readqc
   python <srcdir>/pipeline_readqc.py make full

Requirements
==============

* fastqc
* fastq_screen
* sickle >= 1.33
* cutadapt >= 1.7.1

Code
====

"""

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import shutil
import sqlite3

# import modules from the CGAT code collection
import CGAT.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineReadqc as PipelineReadqc
import CGATPipelines.PipelinePreprocess as PipelinePreprocess
import CGAT.IOTools as IOTools

# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])
PARAMS = P.PARAMS

# define input files and preprocessing steps
# list of acceptable input formats
INPUT_FORMATS = ["*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz"]

# Regular expression to extract a track from an input file. Does not preserve
# a directory as part of the track.
REGEX_TRACK = regex(r"([^/]+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

# Regular expression to extract a track from both processed and unprocessed
# files
REGEX_TRACK_BOTH = regex(r"(processed.dir/)*([^/]+)\.(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(?P<suffix>fastq.1.gz|fastq.gz|sra|csfasta.gz)")

# List of unprocessed input files
UNPROCESSED_INPUT_GLOB = INPUT_FORMATS

# List of processed (trimmed, etc) input files
PROCESSED_INPUT_GLOB = [os.path.join("processed.dir", x)
                        for x in INPUT_FORMATS]


def connect():
    '''
    Setup a connection to an sqlite database
    '''

    dbh = sqlite3.connect(PARAMS['database'])
    return dbh

# if preprocess tools are specified, preprocessing is done on output that has
# already been generated in the first run
if PARAMS.get("preprocessors", None):
    PREPROCESSTOOLS = [tool for tool
                       in P.asList(PARAMS["preprocessors"])]
    preprocess_prefix = ("-".join(PREPROCESSTOOLS[::-1]) + "-")
    if PARAMS["auto_remove"]:
        @follows(mkdir("fasta.dir"))
        @transform(INPUT_FORMATS,
                   SEQUENCEFILES_REGEX,
                   r"fasta.dir/\1.fasta")
        def makeAdaptorFasta(infile, outfile):

            '''
            Make a single fasta file for each sample of all contaminant adaptor
            sequences for removal
            '''
            contams = PARAMS['contaminants']
            PipelinePreprocess.makeAdaptorFasta(infile=infile,
                                                dbh=connect(),
                                                contaminants_file=contams,
                                                outfile=outfile)

        @follows(makeAdaptorFasta)
        @collate(makeAdaptorFasta,
                 regex("fasta.dir/(.+).fasta"),
                 r"%s" % PARAMS['adapter_file'])
        def aggregateAdaptors(infiles, outfile):
            '''
            Collate fasta files into a single contaminants file for
            adapter removal.
            '''

            PipelinePreprocess.mergeAdaptorFasta(infiles, outfile)

    else:
        @follows(mkdir("fasta.dir"))
        @transform(INPUT_FORMATS,
                   SEQUENCEFILES_REGEX,
                   r"fasta.dir/\1.fasta")
        def aggregateAdaptors(infile, outfile):

            P.touch(outfile)

    @follows(mkdir("processed.dir"),
             mkdir("log.dir"),
             mkdir("summary.dir"),
             aggregateAdaptors)
    @transform(INPUT_FORMATS,
               SEQUENCEFILES_REGEX,
               r"processed.dir/%s\1.\g<suffix>" % preprocess_prefix)
    def processReads(infile, outfile):
        '''process reads from .fastq format files
        Tasks specified in PREPROCESSTOOLS are run in order
        '''
        trimmomatic_options = PARAMS["trimmomatic_options"]
        if PARAMS["adapter_file"] or PARAMS["trimmomatic_adapter"]:
            if PARAMS["adapter_file"]:
                adapter_file = PARAMS["adapter_file"]
            else:
                adapter_file = PARAMS["trimmomatic_adapter"]

            adapter_options = " ILLUMINACLIP:%s:%s:%s:%s " % (
                adapter_file,
                PARAMS["trimmomatic_mismatches"],
                PARAMS["trimmomatic_p_thresh"], PARAMS["trimmomatic_c_thresh"])
            trimmomatic_options = adapter_options + trimmomatic_options

        job_threads = PARAMS["threads"]
        job_memory = "7G"

        m = PipelinePreprocess.MasterProcessor(
            save=PARAMS["save"],
            summarise=PARAMS["summarise"],
            threads=PARAMS["threads"],
            trimgalore_options=PARAMS["trimgalore_options"],
            trimmomatic_options=trimmomatic_options,
            sickle_options=PARAMS["sickle_options"],
            flash_options=PARAMS["flash_options"],
            fastx_trimmer_options=PARAMS["fastx_trimmer_options"],
            cutadapt_options=PARAMS["cutadapt_options"],
            adapter_file=PARAMS['adapter_file'])

        statement = m.build((infile,), outfile, PREPROCESSTOOLS)

        P.run()

else:
    @follows(mkdir("processed.dir"))
    def processReads():
        pass


@follows(processReads, mkdir(PARAMS["exportdir"]),
         mkdir(os.path.join(PARAMS["exportdir"], "fastqc")))
@transform(UNPROCESSED_INPUT_GLOB + PROCESSED_INPUT_GLOB,
           REGEX_TRACK,
           r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''run Fastqc on each input file.

    convert sra files to fastq and check mapping qualities are in
    solexa format.  Perform quality control checks on reads from
    .fastq files.
    '''
    # MM: only pass the contaminants file list if requested by user,
    # do not make this the default behaviour
    if PARAMS['add_contaminants']:
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"] + "/fastqc",
                                   contaminants=PARAMS['contaminants'])
    else:
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"] + "/fastqc")

    statement = m.build((infiles,), outfile)
    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(runFastqc, suffix(".fastqc"), "_fastqc.load")
def loadFastqc(infile, outfile):
    '''load FASTQC stats into database.'''
    track = P.snip(infile, ".fastqc")
    filename = os.path.join(
        PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")

    PipelineReadqc.loadFastqc(filename,
                              backend=PARAMS["database_backend"],
                              database=PARAMS["database_name"],
                              host=PARAMS["database_host"],
                              username=PARAMS["database_username"],
                              password=PARAMS["database_password"],
                              port=PARAMS["database_port"])
    P.touch(outfile)


@follows(processReads, mkdir(PARAMS["exportdir"]),
         mkdir(os.path.join(PARAMS["exportdir"], "fastq_screen")))
@transform(UNPROCESSED_INPUT_GLOB + PROCESSED_INPUT_GLOB,
           REGEX_TRACK_BOTH,
           r"%s/fastq_screen/\2.fastq.1_screen.png" % PARAMS['exportdir'])
def runFastqScreen(infiles, outfile):
    '''run FastqScreen on input files.'''

    # variables required for statement built by FastqScreen()
    tempdir = P.getTempDir(".")
    outdir = os.path.join(PARAMS["exportdir"], "fastq_screen")

    # Create fastq_screen config file in temp directory
    # using parameters from Pipeline.ini
    with IOTools.openFile(os.path.join(tempdir, "fastq_screen.conf"),
                          "w") as f:
        for i, k in PARAMS.items():
            if i.startswith("fastq_screen_database"):
                f.write("DATABASE\t%s\t%s\n" % (i[22:], k))

    m = PipelineMapping.FastqScreen()
    statement = m.build((infiles,), outfile)
    P.run()

    shutil.rmtree(tempdir)


@merge(runFastqc, "status_summary.tsv.gz")
def buildFastQCSummaryStatus(infiles, outfile):
    '''load fastqc status summaries into a single table.'''
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildFastQCSummaryStatus(infiles, outfile, exportdir)


@follows(loadFastqc)
@merge(runFastqc, "basic_statistics_summary.tsv.gz")
def buildFastQCSummaryBasicStatistics(infiles, outfile):
    '''load fastqc summaries into a single table.'''
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildFastQCSummaryBasicStatistics(infiles, outfile,
                                                     exportdir)


@follows(mkdir("experiment.dir"), loadFastqc)
@collate(runFastqc,
         regex("(processed.dir/)*(.*)-([^-]*).fastqc"),
         r"experiment.dir/\2_per_sequence_quality.tsv")
def buildExperimentLevelReadQuality(infiles, outfile):
    """
    Collate per sequence read qualities for all replicates per experiment.
    Replicates are the last part of a filename, eg. Experiment-R1,
    Experiment-R2, etc.
    """
    exportdir = os.path.join(PARAMS["exportdir"], "fastqc")
    PipelineReadqc.buildExperimentReadQuality(infiles, outfile, exportdir)


@collate(buildExperimentLevelReadQuality,
         regex("(.+)/(.+)_per_sequence_quality.tsv"),
         r"\1/experiment_per_sequence_quality.tsv")
def combineExperimentLevelReadQualities(infiles, outfile):
    """
    Combine summaries of read quality for different experiments
    """
    infiles = " ".join(infiles)
    statement = ("python %(scriptsdir)s/combine_tables.py "
                 "  --log=%(outfile)s.log "
                 "  --regex-filename='.+/(.+)_per_sequence_quality.tsv' "
                 "%(infiles)s"
                 "> %(outfile)s")
    P.run()


@transform(combineExperimentLevelReadQualities,
           regex(".+/(.+).tsv"),
           r"\1.load")
def loadExperimentLevelReadQualities(infile, outfile):
    P.load(infile, outfile)


@transform((buildFastQCSummaryStatus, buildFastQCSummaryBasicStatistics),
           suffix(".tsv.gz"), ".load")
def loadFastqcSummary(infile, outfile):
    P.load(infile, outfile, options="--add-index=track")


@follows(loadFastqc,
         loadFastqcSummary,
         loadExperimentLevelReadQualities,
         runFastqScreen)
def full():
    pass


@follows(processReads)
def test():
    pass


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


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

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
