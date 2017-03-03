"""
====================
RNASeqQC pipeline
====================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========

This pipeline should be run as the first step in your RNA seq analysis
work flow. It will help detect error and biases within your raw
data. The output of the pipeline can be used to filter out problematic
cells in a standard RNA seq experiment. For single cell RNA seq the
pipeline_rnaseqqc.py should be run instead.

Sailfish is used to perform rapid alignment-free transcript
quantification and hisat is used to align a subset of reads to the
reference genome.

From the sailfish and hisat output, a number of analyses are
performed, either within the pipeline or during the reporting:

- Proportion of reads aligned to annotated features
    (rRNA, protein coding, lincRNA etc)
- Sequencing depth saturation curves Per Sample
- Per-sample expression distributions
- Strandedness assesment
- Assessment of sequence biases
- Expression of top genes and expression of genes of interest

Most of the above analysis will group samples by the sample factors
(see Important configuration options below for details on how factors
are identified)


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.


Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The following suffixes/file types are possible:

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

Important configuration options
===============================

To determine the experimental factors in your experiment, name files
with factors separated by ``-``, for example::

   sample1-mRNA-10k-R1-L01.fastq.1.gz
   sample1-mRNA-10k-R1-L01.fastq.2.gz
   sample1-mRNA-10k-R1-L02.fastq.1.gz
   sample1-mRNA-10k-R1-L02.fastq.2.gz
   sample1-mRNA-150k-R1-L01.fastq.1.gz
   sample1-mRNA-150k-R1-L01.fastq.2.gz
   sample1-mRNA-150k-R1-L02.fastq.1.gz
   sample1-mRNA-150k-R1-L02.fastq.2.gz

and then set the ``factors`` variable in :file:`pipeline.ini` to::

   factors=experiment-source-replicate-lane

If you want to include additional factors which are not identifiable
from the sample names you can specfify these in an optional file
"additional_factors.tsv". This file must contain the sample names in
the first columns and then an additional column for each factor (tab
separated). See below for an example to include the additional factors
"preparation_date" and "rna_quality":

sample    preparation_date    rna_quality
sample1-mRNA-10k-R1-L01    01-01-2016    poor
sample1-mRNA-10k-R1-L01    01-01-2016    poor
sample1-mRNA-10k-R1-L02    04-01-2016    good
sample1-mRNA-10k-R1-L02    04-01-2016    good


Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the
apparent biases in transcript abudance within the sequence archive
The pipeline also produces output in the database file:`csvdb`.

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqqc.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz
   tar -xvzf pipeline_readqc.tgz
   cd pipeline_readqc
   python <srcdir>/pipeline_readqc.py make full

Requirements:

+---------+------------+------------------------------------------------+
|*Program*|*Version*   |*Purpose*                                       |
+---------+------------+------------------------------------------------+
|sailfish |>=0.9.0     |pseudo alignment                               |
+---------+------------+------------------------------------------------+
|hisat    |>=0.1.6     |read mapping                                   |
+---------+------------+------------------------------------------------+
|samtools |>=0.1.16    |bam/sam files
+---------+------------+------------------------------------------------+
|bedtools |            |work with intervals
+---------+------------+------------------------------------------------+
|picard   |>=1.42      |bam/sam files
+---------+------------+------------------------------------------------+
|bamstats |>=1.22      |from CGR, liverpool
+---------+------------+------------------------------------------------+
|sra-tools|            |extracting sra files
+---------+------------+------------------------------------------------+


Glossary
========

.. glossary::

   hisat
      hisat_- a read mapper used in the pipeline because it is
              relatively quick to run
   sailfish
      sailfish_-a pseudoaligner that is used for quantifying the
                abundance transcripts
.._hisat: http://ccb.jhu.edu/software/hisat/manual.shtml
.. sailfish: https://github.com/kingsfordgroup/sailfish

Code
====

"""

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import transform, suffix, regex, merge, \
    follows, mkdir, originate, add_inputs, jobs_limit, split

# import useful standard python modules
import sys
import os
import sqlite3
import re
import pandas as pd
import numpy as np
import itertools
from scipy.stats import linregress

from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineWindows as PipelineWindows
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineRnaseq as PipelineRnaseq

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))


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


#########################################################################
#########################################################################
#########################################################################
# define input files
SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r"(.*\/)*(\S+).(fastq.1.gz|fastq.gz|fa.gz|sra|"
    "csfasta.gz|csfasta.F3.gz|export.txt.gz)")

###################################################################
# Pipeline Utility functions
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])

    if not os.path.exists(PARAMS["annotations_database"]):
        raise ValueError(
            "can't find database '%s'" %
            PARAMS["annotations_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
                (PARAMS["annotations_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


def findSuffixedFile(prefix, suffixes):
    for check_suffix in suffixes:
        check_infile = prefix + check_suffix
        if os.path.exists(check_infile):
            return (check_infile, check_suffix)

###################################################################
# count number of reads
###################################################################


@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\2.nreads")
def countReads(infile, outfile):
    '''Count number of reads in input files.'''
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()

###################################################################
# build geneset
###################################################################


@follows(mkdir("geneset.dir"))
@merge(PARAMS["annotations_interface_geneset_all_gtf"],
       "geneset.dir/reference.gtf.gz")
def buildReferenceGeneSet(infile, outfile):
    ''' filter full gene set and add attributes to create the reference gene set

    Performs merge and filter operations:
       * Merge exons separated by small introns (< 5bp).
       * Remove transcripts with very long introns (`max_intron_size`)
       * Remove transcripts located on contigs to be ignored (`remove_contigs`)
         (usually: chrM, _random, ...)
       * (Optional) Remove transcripts overlapping repetitive sequences
         (`rna_file`)

    This preserves all features in a gtf file (exon, CDS, ...)

    Runs cuffcompare with `infile` against itself to add
    attributes such as p_id and tss_id.

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile : str
       Input filename in :term:`gtf` format
    annotations_interface_rna_gff : str
       :term:`PARAMS`. Filename of :term:`gtf` file containing
       repetitive rna annotations
    genome_dir : str
       :term:`PARAMS`. Directory of :term:fasta formatted files
    genome : str
       :term:`PARAMS`. Genome name (e.g hg38)
    '''

    tmp_mergedfiltered = P.getTempFilename(".")

    if "geneset_remove_repetetive_rna" in PARAMS:
        rna_file = PARAMS["annotations_interface_rna_gff"]
    else:
        rna_file = None

    gene_ids = PipelineMapping.mergeAndFilterGTF(
        infile,
        tmp_mergedfiltered,
        "%s.removed.gz" % outfile,
        genome=os.path.join(PARAMS["genome_dir"], PARAMS["genome"]),
        max_intron_size=PARAMS["max_intron_size"],
        remove_contigs=PARAMS["geneset_remove_contigs"],
        rna_file=rna_file)

    # Add tss_id and p_id
    PipelineMapping.resetGTFAttributes(
        infile=tmp_mergedfiltered,
        genome=os.path.join(PARAMS["genome_dir"], PARAMS["genome"]),
        gene_ids=gene_ids,
        outfile=outfile)

    os.unlink(tmp_mergedfiltered)


@follows(mkdir("geneset.dir"))
@originate("geneset.dir/protein_coding_gene_ids.tsv")
def identifyProteinCodingGenes(outfile):
    '''Output a list of proteing coding gene identifiers

    Identify protein coding genes from the annotation database table
    and output the gene identifiers

    Parameters
    ----------
    oufile : str
       Output file of :term:`gtf` format
    annotations_interface_table_gene_info : str
       :term:`PARAMS`. Database table name for gene information

    '''

    dbh = connect()

    table = os.path.basename(PARAMS["annotations_interface_table_gene_info"])

    select = dbh.execute("""SELECT DISTINCT gene_id
    FROM annotations.%(table)s
    WHERE gene_biotype = 'protein_coding'""" % locals())

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("gene_id\n")
        outf.write("\n".join((x[0] for x in select)) + "\n")


@transform(buildReferenceGeneSet,
           suffix("reference.gtf.gz"),
           add_inputs(identifyProteinCodingGenes),
           "refcoding.gtf.gz")
def buildCodingGeneSet(infiles, outfile):
    '''build a gene set with only protein coding transcripts.

    Retain the genes from the gene_tsv file in the outfile geneset.
    The gene set will contain all transcripts of protein coding genes,
    including processed transcripts. The gene set includes UTR and
    CDS.

    Parameters
    ----------
    infiles : list
    infile: str
       Input filename in :term:`gtf` format

    genes_ts: str
       Input filename in :term:`tsv` format

    outfile: str
       Output filename in :term:`gtf` format

    '''

    infile, genes_tsv = infiles

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter
    --filter-method=gene
    --map-tsv-file=%(genes_tsv)s
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


@follows(mkdir("geneset.dir"))
@merge(PARAMS["annotations_interface_geneset_all_gtf"],
       "geneset.dir/coding_exons.gtf.gz")
def buildCodingExons(infile, outfile):
    '''compile the set of protein coding exons.

    Filter protein coding transcripts
    This set is used for splice-site validation

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile: str
       Output filename in :term:`gtf` format

    '''

    statement = '''
    zcat %(infile)s
    | awk '$3 == "CDS"'
    | cgat gtf2gtf
    --method=filter
    --filter-method=proteincoding
    --log=%(outfile)s.log
    | awk -v OFS="\\t" -v FS="\\t" '{$3="exon"; print}'
    | cgat gtf2gtf
    --method=merge-exons
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


@transform(buildCodingGeneSet, suffix(".gtf.gz"), ".junctions")
def buildJunctions(infile, outfile):
    '''build file with splice junctions from gtf file.

    Identify the splice junctions from a gene set :term:`gtf`
    file. A junctions file is a better option than supplying a GTF
    file, as parsing the latter often fails. See:

    http://seqanswers.com/forums/showthread.php?t=7563

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile: str
       Output filename

    '''

    outf = IOTools.openFile(outfile, "w")
    njunctions = 0
    for gffs in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile, "r"))):

        gffs.sort(key=lambda x: x.start)
        end = gffs[0].end
        for gff in gffs[1:]:
            # subtract one: these are not open/closed coordinates but
            # the 0-based coordinates
            # of first and last residue that are to be kept (i.e., within the
            # exon).
            outf.write("%s\t%i\t%i\t%s\n" %
                       (gff.contig, end - 1, gff.start, gff.strand))
            end = gff.end
            njunctions += 1

    outf.close()

    if njunctions == 0:
        E.warn('no junctions found in gene set')
        return
    else:
        E.info('found %i junctions before removing duplicates' % njunctions)

    # make unique
    statement = '''mv %(outfile)s %(outfile)s.tmp;
                   cat < %(outfile)s.tmp | sort | uniq > %(outfile)s;
                   rm -f %(outfile)s.tmp; '''
    P.run()


@transform(buildCodingGeneSet,
           suffix(".gtf.gz"),
           ".fasta")
def buildTranscriptFasta(infile, outfile):
    """build geneset where all exons within a gene
    are merged.
    """
    dbname = outfile[:-len(".fasta")]

    statement = '''zcat %(infile)s
    | cgat gff2fasta
    --is-gtf
    --genome=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | cgat index_fasta
    %(dbname)s --force-output -
    > %(dbname)s.log
    '''
    P.run()


@transform(buildCodingGeneSet,
           suffix(".gtf.gz"),
           ".tsv")
def buildTranscriptGeneMap(infile, outfile):
    """build a map of transcript ids to gene ids."""

    statement = """
    zcat %(infile)s
    |cgat gtf2tsv
    --attributes-as-columns
    --output-only-attributes
    | cgat csv_cut transcript_id gene_id
    > %(outfile)s"""
    P.run()

###################################################################
# subset fastqs
###################################################################


@follows(mkdir("fastq.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"fastq.dir/\2.subset")
def subsetSequenceData(infile, outfile):
    """subset fastq files"""
    ignore_pipe_erors = True
    ignore_errors = True
    m = PipelineMapping.SubsetHead(limit=PARAMS["sample_size"])
    statement = m.build((infile,), outfile)
    P.run()
    P.touch(outfile)


@follows(mkdir("fastq.dir"))
@merge(countReads,
       "fastq.dir/highest_depth_sample.sentinel")
def identifyHighestDepth(infiles, outfile):
    ''' identify the sample with the highest depth'''

    highest_depth = 0
    for count_inf in infiles:
        for line in IOTools.openFile(count_inf, "r"):
            if not line.startswith("nreads"):
                continue
            nreads = int(line[:-1].split("\t")[1])
            if nreads > highest_depth:
                highest_depth = nreads
                highest_depth_sample = os.path.basename(
                    P.snip(count_inf, ".nreads"))

    assert highest_depth_sample, ("unable to identify the sample "
                                  "with the highest depth")

    infile, inf_suffix = findSuffixedFile(highest_depth_sample,
                                          [x[1:] for x in SEQUENCESUFFIXES])
    infile = os.path.abspath(infile)

    assert infile, ("unable to find the raw data for the "
                    "sample with the highest depth")

    dst = os.path.abspath(P.snip(outfile, ".sentinel") + inf_suffix)

    def forcesymlink(src, dst):
        try:
            os.symlink(src, dst)
        except:
            os.remove(dst)
            os.symlink(src, dst)

    forcesymlink(infile, dst)

    # if paired end fastq, need to link the paired end too!
    if inf_suffix == ".fastq.1.gz":
        dst2 = P.snip(outfile, ".sentinel") + ".fastq.2.gz"
        forcesymlink(infile.replace(".fastq.1.gz", ".fastq.2.gz"), dst2)

    forcesymlink("%s.nreads" % highest_depth_sample,
                 "nreads.dir/highest_depth_sample.nreads")

    P.touch(outfile)


@split(identifyHighestDepth,
       "fastq.dir/highest_counts_subset_*")
def subsetRange(infile, outfiles):
    '''subset highest depth sample to 10%-100% depth '''

    outfile = "fastq.dir/highest_counts_subset.sentinel"
    infile_prefix = P.snip(os.path.basename(infile), ".sentinel")
    nreads_inf = "nreads.dir/%s.nreads" % infile_prefix

    for line in IOTools.openFile(nreads_inf, "r"):
        if not line.startswith("nreads"):
            continue
        nreads = int(line[:-1].split("\t")[1])

    infile, inf_suffix = findSuffixedFile(P.snip(infile, ".sentinel"),
                                          [x[1:] for x in SEQUENCESUFFIXES])

    # PipelineMapping.Counter double counts for paired end
    # Note: this wont handle sra. Need to add a call to Sra.peak to check for
    # paired end files in SRA
    if inf_suffix == ".fastq.1.gz":
        nreads = nreads / 2

    subset_depths = list(range(10, 110, 10))
    limits = [int(nreads / (100.0 / int(depth)))
              for depth in subset_depths]

    ignore_pipe_erors = True
    ignore_errors = True
    m = PipelineMapping.SubsetHeads(limits=limits)
    statement = m.build((infile,), outfile)

    P.run()

    P.touch(outfile)


@follows(subsetSequenceData)
def subset():
    pass

###################################################################
# map reads
###################################################################


@follows(mkdir("hisat.dir"))
@transform(subsetSequenceData,
           regex("fastq.dir/(.*).subset"),
           add_inputs(buildJunctions),
           r"hisat.dir/\1.hisat.bam")
def mapReadsWithHisat(infiles, outfile):
    '''
    Map reads using Hisat  (spliced reads).

    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        filename with suffix .junctions containing a list of known
        splice junctions.

    hisat_threads: int
        :term:`PARAMS`
        number of threads with which to run hisat

    hisat_memory: str
        :term:`PARAMS`
        memory required for hisat job

    hisat_executable: str
        :term:`PARAMS`
        path to hisat executable

    hisat_library_type: str
        :term:`PARAMS`
        hisat rna-strandess parameter, see
        https://ccb.jhu.edu/software/hisat/manual.shtml#command-line

    hisat_options: str
        options string for hisat, see
        https://ccb.jhu.edu/software/hisat/manual.shtml#command-line

    hisat_index_dir: str
        path to directory containing hisat indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.

    .. note::
    If hisat fails with an error such as::

       Error: segment-based junction search failed with err =-6
       what():  std::bad_alloc

    it means that it ran out of memory.

    '''

    job_threads = PARAMS["hisat_threads"]
    job_memory = PARAMS["hisat_memory"]

    m = PipelineMapping.Hisat(
        executable=P.substituteParameters(
            **locals())["hisat_executable"],
        strip_sequence=PARAMS["strip_sequence"])

    infile, junctions = infiles
    infile = P.snip(infile, ".subset") + ".fastq.gz"
    if not os.path.exists(infile):
        infile = P.snip(infile, ".fastq.gz") + ".fastq.1.gz"

    statement = m.build((infile,), outfile)

    P.run()


###################################################################
# build mapping stats
###################################################################


@transform(mapReadsWithHisat,
           regex("(.*)/(.*)\.(.*).bam"),
           r"\1/\2.\3.readstats")
def buildBAMStats(infile, outfile):
    '''count number of reads mapped, duplicates, etc.

    Excludes regions overlapping repetitive RNA sequences

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename with number of reads per sample

    outfile : str
       Output filename with read stats

    annotations_interface_rna_gtf : str
        :term:`PARMS`. :term:`gtf` format file with repetitive rna
    '''

    rna_file = PARAMS["annotations_interface_rna_gff"]

    job_memory = "16G"

    track = P.snip(os.path.basename(infile), ".hisat.bam")

    # if a fastq file exists, submit for counting
    if os.path.exists(track + ".fastq.gz"):
        fastqfile = track + ".fastq.gz"
    elif os.path.exists(track + ".fastq.1.gz"):
        fastqfile = track + ".fastq.1.gz"
    else:
        fastqfile = None

    if fastqfile is not None:
        fastq_option = "--fastq-file=%s" % fastqfile
    else:
        fastq_option = ""

    statement = '''
    cgat bam2stats
         %(fastq_option)s
         --force-output
         --mask-bed-file=%(rna_file)s
         --ignore-masked-reads
         --num-reads=%(sample_size)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()


@P.add_doc(PipelineMappingQC.loadBAMStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildBAMStats, "bam_stats.load")
def loadBAMStats(infiles, outfile):
    ''' load bam statistics into bam_stats table '''
    PipelineMappingQC.loadBAMStats(infiles, outfile)


@P.add_doc(PipelineWindows.summarizeTagsWithinContext)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(
               PARAMS["annotations_interface_genomic_context_bed"]),
           ".contextstats.tsv.gz")
def buildContextStats(infiles, outfile):
    ''' build mapping context stats '''
    PipelineWindows.summarizeTagsWithinContext(
        infiles[0], infiles[1], outfile)


@P.add_doc(PipelineWindows.loadSummarizedContextStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadBAMStats)
@merge(buildContextStats, "context_stats.load")
def loadContextStats(infiles, outfile):
    ''' load context mapping statistics into context_stats table '''
    PipelineWindows.loadSummarizedContextStats(infiles, outfile)


@originate("geneset.dir/altcontext.bed.gz")
def buildBedContext(outfile):
    ''' Generate a bed file that can be passed into buildAltContextStats '''

    dbh = connect()

    tmp_bed_sorted_filename = P.getTempFilename(shared=True)

    tmp_bed_sorted = IOTools.openFile(tmp_bed_sorted_filename, "w")

    sql_statements = [
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "lincRNA"
        FROM gene_info GI
        JOIN geneset_lincrna_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "lincRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "snoRNA"
        FROM gene_info GI
        JOIN geneset_noncoding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "snoRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "miRNA"
        FROM gene_info GI
        JOIN geneset_noncoding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "miRNA"''',
        '''SELECT DISTINCT GTF.contig, GTF.start, GTF.end, "protein_coding"
        FROM gene_info GI
        JOIN geneset_coding_exons_gtf GTF
        ON GI.gene_id=GTF.gene_id
        WHERE GI.gene_biotype == "protein_coding"''']

    for sql_statement in sql_statements:
        state = dbh.execute(sql_statement)

        for line in state:
            tmp_bed_sorted.write(("%s\n") % "\t".join(map(str, line)))

    tmp_bed_sorted.close()

    statement = '''sortBed  -i %(tmp_bed_sorted_filename)s
    | gzip > %(outfile)s'''

    P.run()

    os.unlink(tmp_bed_sorted_filename)


@P.add_doc(PipelineWindows.summarizeTagsWithinContext)
@follows(buildBedContext)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(buildBedContext),
           ".altcontextstats.tsv.gz")
def buildAltContextStats(infiles, outfile):
    ''' build mapping context stats of snoRNA, miRNA,
        lincRNA, protein coding '''

    infile, bed = infiles

    PipelineWindows.summarizeTagsWithinContext(
        infile, bed,  outfile)


@P.add_doc(PipelineWindows.loadSummarizedContextStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(loadContextStats)
@merge(buildAltContextStats, "altcontext_stats.load")
def loadAltContextStats(infiles, outfile):
    ''' load context mapping statistics into context_stats table '''
    PipelineWindows.loadSummarizedContextStats(infiles,
                                               outfile,
                                               suffix=".altcontextstats.tsv.gz")


###################################################################
# alignment-free quantification
###################################################################


@follows(mkdir("sailfish.dir"))
@transform(buildTranscriptFasta,
           regex("(\S+)"),
           "sailfish.dir/transcripts.sailfish.index")
def indexForSailfish(infile, outfile):
    '''create a sailfish index'''

    statement = '''
    sailfish index --transcripts=%(infile)s
    --out=%(outfile)s '''
    P.run()


@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(indexForSailfish,
                      buildCodingGeneSet,
                      buildTranscriptGeneMap),
           r"sailfish.dir/\2/quant.sf")
def runSailfish(infiles, outfile):
    '''quantify abundance'''

    job_threads = PARAMS["sailfish_threads"]
    job_memory = PARAMS["sailfish_memory"]

    infile, index, geneset, transcript_map = infiles

    sailfish_bootstrap = 1
    sailfish_libtype = PARAMS["sailfish_libtype"]
    sailfish_options = PARAMS["sailfish_options"]
    sailfish_options += " --geneMap %s" % transcript_map

    m = PipelineMapping.Sailfish()

    statement = m.build((infile,), outfile)

    P.run()


@split(runSailfish,
       ["sailfish.dir/sailfish_transcripts.tsv.gz",
        "sailfish.dir/sailfish_genes.tsv.gz"])
def mergeSailfishResults(infiles, outfiles):
    ''' concatenate sailfish expression estimates from each sample'''

    s_infiles = " " .join(sorted(infiles))
    outfile_transcripts, outfile_genes = outfiles

    statement = """
    cat %(s_infiles)s
    | awk -v OFS="\\t"
    '/^Name/
    { sample_id+=1;
      if (sample_id == 1) {
         gsub(/Name/, "transcript_id");
         printf("sample_id\\t%%s\\n", $0)};
      next;}
    !/^#/
        {printf("%%i\\t%%s\\n", sample_id, $0)}'
    | gzip
    > %(outfile_transcripts)s
    """
    P.run()

    s_infiles = " ".join(
        [re.sub("quant.sf", "quant.genes.sf", x) for x in infiles])

    statement = """
    cat %(s_infiles)s
    | awk -v OFS="\\t"
    '/^Name/
    { sample_id+=1;
      if (sample_id == 1) {
         gsub(/Name/, "gene_id");
         printf("sample_id\\t%%s\\n", $0)};
      next;}
    !/^#/
        {printf("%%i\\t%%s\\n", sample_id, $0)}'
    | gzip
    > %(outfile_genes)s
    """
    P.run()


@transform(mergeSailfishResults,
           suffix(".tsv.gz"),
           ".load")
def loadSailfishResults(infile, outfile):
    P.load(infile, outfile,
           options="--add-index=sample_id "
           "--add-index=gene_id "
           "--add-index=transcript_id "
           "--map=sample_id:int")

###################################################################
# strand bias
###################################################################


@transform(buildReferenceGeneSet,
           suffix("reference.gtf.gz"),
           "refflat.txt")
def buildRefFlat(infile, outfile):
    '''build flat geneset for Picard RnaSeqMetrics.'''

    tmpflat = P.getTempFilename(".")

    statement = '''
    gtfToGenePred -genePredExt -geneNameAsName2 %(infile)s %(tmpflat)s;
    paste <(cut -f 12 %(tmpflat)s) <(cut -f 1-10 %(tmpflat)s)
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmpflat)


@P.add_doc(PipelineMappingQC.buildPicardRnaSeqMetrics)
@transform(mapReadsWithHisat,
           suffix(".bam"),
           add_inputs(buildRefFlat),
           ".picard_rna_metrics")
def buildPicardRnaSeqMetrics(infiles, outfile):
    '''Get duplicate stats from picard RNASeqMetrics '''
    # convert strandness to tophat-style library type
    if PARAMS["hisat_library_type"] == ("RF" or "R"):
        strand = "SECOND_READ_TRANSCRIPTION_STRAND"
    elif PARAMS["hisat_library_type"] == ("FR" or "F"):
        strand = "FIRST_READ_TRANSCRIPTION_STRAND"
    else:
        strand = "NONE"

    PipelineMappingQC.buildPicardRnaSeqMetrics(infiles, strand, outfile)


@P.add_doc(PipelineMappingQC.loadPicardRnaSeqMetrics)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardRnaSeqMetrics, ["picard_rna_metrics.load",
                                  "picard_rna_histogram.load"])
def loadPicardRnaSeqMetrics(infiles, outfiles):
    '''merge alignment stats into single tables.'''
    PipelineMappingQC.loadPicardRnaSeqMetrics(infiles, outfiles)


###################################################################
# saturation analysis
###################################################################

@transform(subsetRange,
           regex("fastq.dir/highest_counts_subset_(\d+)."
                 "(fastq.1.gz|fastq.gz|fa.gz|sra|"
                 "csfasta.gz|csfasta.F3.gz|export.txt.gz)"),
           add_inputs(indexForSailfish,
                      buildCodingGeneSet,
                      buildTranscriptGeneMap),
           r"sailfish.dir/highest_counts_subset_\1/quant.sf")
def runSailfishSaturation(infiles, outfile):
    '''quantify abundance of transcripts with increasing subsets of the data'''

    job_threads = PARAMS["sailfish_threads"]
    job_memory = PARAMS["sailfish_memory"]

    infile, index, geneset, transcript_map = infiles

    sailfish_bootstrap = 20
    sailfish_libtype = PARAMS["sailfish_libtype"]
    sailfish_options = PARAMS["sailfish_options"]
    sailfish_options += " --geneMap %s" % transcript_map

    m = PipelineMapping.Sailfish()

    statement = m.build((infile,), outfile)

    P.run()


@mkdir("sailfish.dir/plots.dir")
@merge(runSailfishSaturation,
       "sailfish.dir/plots.dir/saturation_plots.sentinel")
def plotSailfishSaturation(infiles, outfile):
    ''' Plot the relationship between sample sequencing depth and
    quantification accuracy'''

    plotfile_base = P.snip(outfile, ".sentinel")
    bootstrap_sat_plotfile = plotfile_base + "_boostrap_cv.png"
    accuracy_sat_plotfile = plotfile_base + "_accuracy.png"

    quant_dir = os.path.dirname(os.path.dirname(infiles[0]))

    # This is currently hardcoded to expect 10 infiles named:
    # (quant.dir)/highest_counts_subset_(n)/quant.sf,
    # where (n) is the subset index (0-9)

    R('''
    library(reshape2)
    library(ggplot2)
    library(Hmisc)

    Path = "%(quant_dir)s"

    # following code to read Sailfish binary files borrows from Rob
    # Patro's Wasabi R package for making sailfish/salmon output
    # compatable with sleuth
    minfo <- rjson::fromJSON(file=file.path(
      Path, 'highest_counts_subset_9', "aux", "meta_info.json"))

    numBoot <- minfo$num_bootstraps

    point_df = read.table(file.path(Path, 'highest_counts_subset_9', "quant.sf"),
                          sep="\t", header=T, row.names=1)

    final_cols = NULL

    for (ix in seq(0,9,1)){
      bootCon <- gzcon(file(file.path(
        Path, paste0('highest_counts_subset_', ix), 'aux',
                     'bootstrap', 'bootstraps.gz'), "rb"))

      # read in binary data
      boots <- readBin(bootCon, "double",
                       n = minfo$num_targets * minfo$num_bootstraps)
      close(bootCon)

      # turn data into dataframe
      boots_df = t(data.frame(matrix(unlist(boots),
                              nrow=minfo$num_bootstraps, byrow=T)))

      # add rownames
      rownames(boots_df) = rownames(point_df)

      final_cols[[paste0("sample_", ix)]] = apply(boots_df, 1,
                                                  function(x) sd(x)/mean(x))
    }

    # make final dataframe with boostrap CVs
    final_df = data.frame(do.call("cbind", final_cols))

    # add expression values, subset to transcripts with >1 read and bin exp
    final_df$max_exp = point_df$NumReads
    final_df = final_df[final_df$max_exp>1,]
    final_df$max_exp = as.numeric(cut2(final_df$max_exp, g=10))

    # melt and aggregate
    melted_df = melt(final_df, id="max_exp")
    melted_df = melted_df[is.finite(melted_df$value),]
    aggdata <-aggregate(melted_df$value,
                        by=list(melted_df$max_exp, melted_df$variable),
                        FUN=mean)
    aggdata$Group.1 = as.factor(aggdata$Group.1)

    m_txt = element_text(size=20)
    my_theme = theme(
    axis.text=m_txt,
    axis.title=m_txt,
    legend.text=m_txt,
    legend.title=m_txt,
    aspect.ratio=1)

    p = ggplot(aggdata, aes(10*as.numeric(Group.2), x,
                            colour=Group.1, group=Group.1)) +
    geom_line() +
    theme_bw() +
    xlab("Sampling depth (%%)") +
    ylab("Average Coefficient of variance") +
    scale_colour_manual(name="Exp. Decile",
                        values=colorRampPalette(c("yellow","purple"))(10)) +
    scale_x_continuous(breaks=seq(10,100,10), limits=c(10,100)) +
    my_theme

    ggsave("%(bootstrap_sat_plotfile)s")

    # read in the point estimate data

    tpm_est = NULL

    ref_point_df = read.table(
      file.path(Path, 'highest_counts_subset_9', "quant.sf"),
      sep="\t", header=T, row.names=1)

    for (ix in seq(0,9,1)){

    point_df = read.table(
      file.path(Path, paste0('highest_counts_subset_', ix), "quant.sf"),
      sep="\t", header=T, row.names=1)

    tpm_est[[paste0("sample_", ix)]] = (
      abs(point_df$TPM - ref_point_df$TPM) / ref_point_df$TPM)
    }

    tpm_est_df = data.frame(do.call("cbind", tpm_est))

    # add expression values, subset to transcripts with >1 read and bin exp.
    tpm_est_df$max_exp = point_df$NumReads
    tpm_est_df = tpm_est_df[point_df$NumReads>1,]
    tpm_est_df$max_exp = as.numeric(cut2(tpm_est_df$max_exp, g=10))

    # melt and aggregate
    melted_df = melt(tpm_est_df, id="max_exp")
    melted_df = melted_df[is.finite(melted_df$value),]
    aggdata <-aggregate(melted_df$value,
                        by=list(melted_df$max_exp, melted_df$variable),
                        FUN=mean)
    aggdata$Group.1 = as.factor(aggdata$Group.1)

    p = ggplot(aggdata, aes(10*as.numeric(Group.2), x,
                            colour=Group.1, group=Group.1)) +
    geom_line() +
    theme_bw() +
    xlab("Sampling depth (%%)") +
    ylab("Abs. difference in exp. estimate (normalised)") +
    scale_colour_manual(name="Exp. Decile",
                        values=colorRampPalette(c("yellow","purple"))(10)) +
    scale_x_continuous(breaks=seq(10,90,10), limits=c(10,90)) +
    my_theme

    ggsave("%(accuracy_sat_plotfile)s")

    ''' % locals())

    P.touch(outfile)


###################################################################
# gene coverage profiles
###################################################################


@follows(mkdir("transcriptprofiles.dir"))
@transform(mapReadsWithHisat,
           regex("hisat.dir/(\S+).hisat.bam"),
           add_inputs(buildCodingExons),
           r"transcriptprofiles.dir/\1.transcriptprofile.gz")
def buildTranscriptProfiles(infiles, outfile):
    '''build gene coverage profiles

    PolyA-RNA-Seq is expected to show a bias towards the 3' end of
    transcripts. Here we generate a meta-profile for each sample for
    the read depth from the :term:`bam` file across the gene models
    defined in the :term:`gtf` gene set

    In addition to the outfile specified by the task, plots will be
    saved with full and focus views of the meta-profile

    Parameters
    ----------
    infiles : list of str
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str`
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`tsv` format
    '''

    bamfile, gtffile = infiles

    job_memory = "8G"

    statement = '''cgat bam2geneprofile
    --output-filename-pattern="%(outfile)s.%%s"
    --force-output
    --reporter=transcript
    --use-base-accuracy
    --method=geneprofile
    --method=geneprofileabsolutedistancefromthreeprimeend
    --normalize-profile=all
    %(bamfile)s %(gtffile)s
    | gzip
    > %(outfile)s
    '''

    P.run()


@merge(buildTranscriptProfiles,
       "transcriptprofiles.dir/threeprimebiasprofiles.load")
def loadTranscriptProfiles(infiles, outfile):
    ''' concatenate and load the transcript profiles
    Retain sample name as column = "track'''

    regex = ("transcriptprofiles.dir/(\S+).transcriptprofile.gz."
             "geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz")

    infiles = [
        x + ".geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz" for x in infiles]

    P.concatenateAndLoad(infiles, outfile, regex_filename=regex)


@merge(SEQUENCEFILES,
       "experiment.tsv")
def buildExperimentTable(infiles, outfile):

    d = os.getcwd()
    try:
        project_id = P.getProjectId()
    except ValueError:
        project_id = "unknown"
    with IOTools.openFile(outfile, "w") as outf:
        outf.write("id\tname\tproject_id\tdirectory\ttitle\n")
        outf.write("\t".join(
            ("1",
             P.getProjectName(),
             project_id,
             d,
             PARAMS.get("title", ""))) + "\n")


@merge(SEQUENCEFILES,
       "samples.tsv")
def buildSamplesTable(infiles, outfile):

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("id\texperiment_id\tsample_name\n")

        for sample_id, filename in enumerate(sorted(infiles)):
            sample_name, suffix = os.path.basename(filename).split(".", 1)
            outf.write("\t".join(
                (str(sample_id + 1), "1", sample_name)) + "\n")


@merge(SEQUENCEFILES,
       "factors.tsv")
def buildFactorTable(infiles, outfile):

    if "factors" not in PARAMS:
        raise ValueError("factors not defined in config file")

    factor_names = PARAMS.get("factors")
    if factor_names is None or factor_names == "!?":
        raise ValueError("factors not defined in config file")
    factor_names = factor_names.split("-")

    sampleID2sampleName = {}

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("sample_id\tfactor\tfactor_value\n")

        for sample_id, filename in enumerate(sorted(infiles)):

            sample_name, suffix = os.path.basename(filename).split(".", 1)
            sampleID2sampleName[sample_name] = sample_id + 1

            parts = sample_name.split("-")

            if len(parts) != len(factor_names):
                raise ValueError(
                    "unexpected number of factors in sample {}: "
                    "expected={}, got={}".format(
                        filename, factor_names, parts))

            for factor, factor_value in zip(factor_names, parts):
                if factor == "_":
                    continue
                outf.write("\t".join((str(sample_id + 1),
                                      factor, factor_value)) + "\n")
            outf.write("\t".join((str(sample_id + 1), "genome",
                                  PARAMS["genome"])) + "\n")

        if os.path.exists("additional_factors.tsv"):
            with IOTools.openFile("additional_factors.tsv", "r") as inf:
                header = next(inf)
                header = header.strip().split("\t")
                additional_factors = header[1:]
                for line in inf:
                    line = line.strip().split("\t")
                    sample_name = line[0]
                    factors_values = line[1:]
                    for factor_ix in range(0, len(additional_factors)):
                        try:
                            outf.write("\t".join((
                                str(sampleID2sampleName[sample_name]),
                                additional_factors[factor_ix],
                                factors_values[factor_ix])) + "\n")
                        except KeyError as ke:
                            sample_names = [os.path.basename(x).split(".")[0]
                                            for x in infiles]
                            raise KeyError(
                                "Sample name in additional_factors table does "
                                " not match up with sample names from raw "
                                "infiles: %s not in %s" % (
                                    ke, ",".join(sample_names)))


@transform((buildExperimentTable, buildSamplesTable, buildFactorTable),
           suffix(".tsv"),
           ".load")
def loadMetaInformation(infile, outfile):
    P.load(infile, outfile,
           options="--map=id:int "
           "--map=sample_id:int "
           "--map=experiment_id:int "
           "--add-index=id "
           "--add-index=experiment_id "
           "--add-index=sample_id ")


@transform(buildTranscriptFasta,
           suffix("refcoding.fasta"),
           "transcripts_attributes.tsv.gz")
def characteriseTranscripts(infile, outfile):
    ''' obtain attributes for transcripts '''

    statement = '''
    cat %(infile)s | cgat fasta2table
    --split-fasta-identifier --section=na,dn,length -v 0
    | gzip > %(outfile)s'''

    P.run()


@transform(characteriseTranscripts,
           regex("transcripts_attributes.tsv.gz"),
           add_inputs(mergeSailfishResults),
           "bias_binned_means.tsv")
def summariseBias(infiles, outfile):

    def percentage(x):
        return float(x[0]) / float(x[1])

    def lin_reg_grad(x, y):
        slope, intercept, r, p, stderr = linregress(x, y)
        return slope

    attributes, genes, transcripts = infiles

    atr = pd.read_table(attributes, sep='\t', index_col="id")
    atr = atr.rename(columns={'pGC': 'GC_Content'})

    for di in itertools.product("ATCG", repeat=2):
        di = di[0] + di[1]
        temp_df = atr.loc[:, [di, "length"]]
        atr[di] = temp_df.apply(percentage, axis=1)

    drop_cols = (["nAT", "nGC", "pAT", "pA", "pG", "pC", "pT", "nA",
                  "nG", "nC", "nT", "ncodons",
                  "mCountsOthers", "nUnk", "nN", "pN"])
    atr = atr.drop(drop_cols, axis=1)

    atr["length"] = np.log2(atr["length"])

    E.info("loading transcripts from {}".format(transcripts))
    exp = pd.read_csv(transcripts, sep='\t', index_col="transcript_id")
    exp['LogTPM'] = np.log2(exp['TPM'] + 0.1)

    merged = atr.join(exp[['sample_id', 'LogTPM']])

    def norm(array):
        array_min = array.min()
        array_max = array.max()
        return pd.Series([(x - array_min) / (array_max - array_min) for x in array])

    def bin2floats(qcut_bin):
        qcut_bin2 = qcut_bin.replace("(", "[").replace(")", "]")
        try:
            qcut_list = eval(qcut_bin2, {'__builtins__': None}, {})
            return qcut_list
        except:
            print("FAILED!!! qcut_bin: ", qcut_bin2)
            return None

    def aggregate_by_factor(df, attribute, sample_names, bins, function):

        temp_dict = dict.fromkeys(sample_names, function)

        temp_dict[attribute] = function
        means_df = df[["LogTPM", "sample_id"]].groupby(
            ["sample_id", pd.qcut(df.ix[:, attribute], bins)])

        means_df = pd.DataFrame(means_df.agg(function))
        means_df.reset_index(inplace=True)

        atr_values = means_df[attribute]
        means_df.drop(attribute, axis=1, inplace=True)

        means_df["LogTPM_norm"] = list(
            means_df.groupby("sample_id")["LogTPM"].apply(norm))

        means_df[attribute] = [np.mean(bin2floats(x)) for x in atr_values]
        means_df = pd.melt(means_df, id_vars=[attribute, "sample_id"])
        means_df.columns = ["bin", "sample_id", "variable", "value"]
        means_df["bias_factor"] = [attribute, ] * len(means_df)

        return means_df

    means_binned_df = pd.DataFrame()

    samples = set(exp.index)
    factors = atr.columns.tolist()

    for factor in factors:
        tmp_df = aggregate_by_factor(
            merged, factor, samples, PARAMS["bias_bin"], np.mean)

        means_binned_df = pd.concat([means_binned_df, tmp_df], axis=0)

    means_binned_df.to_csv(outfile, sep="\t",
                           index=False, float_format='%.6f')


@transform(summariseBias,
           suffix(".tsv"),
           ".load")
def loadBias(infile, outfile):
    P.load(infile, outfile, options="--add-index=sample_id")


###################################################################
# top genes
###################################################################


@mkdir("sailfish.dir/plots.dir")
@follows(loadSailfishResults, loadMetaInformation)
@originate("sailfish.dir/plots.dir/top_expressed.sentinel")
def plotTopGenesHeatmap(outfile):
    '''extract the top 1000 genes (by expression) for each sample and
    plot a heatmap of the intersection'''

    # if someone can find a nice heatmap plotter from a dissimilarity
    # matrix which is compatable with CGATReport, the sqlite and
    # pandas code should be changed into a tracker

    exp_select_cmd = '''
    SELECT TPM, gene_id, sample_name
    FROM sailfish_genes AS A
    JOIN samples AS B
    ON A.sample_id = B.id
    '''

    dbh = connect()

    exp_df = pd.read_sql(exp_select_cmd, dbh)

    factors_select_cmd = '''
    SELECT factor, factor_value, sample_name
    FROM samples AS A
    JOIN factors AS B
    ON A.id = B.sample_id
    '''

    top_n = 1000

    factors_df = pd.read_sql(factors_select_cmd, dbh)

    exp_df['TPM'] = exp_df['TPM'].astype(float)
    exp_df_pivot = pd.pivot_table(exp_df, values=["TPM"],
                                  index="gene_id",
                                  columns="sample_name")

    # extract the top genes per sample
    top_genes = {}
    for col in exp_df_pivot.columns:
        top_genes[col] = exp_df_pivot[col].sort_values(
            ascending=False)[0:top_n].index

    # set up the empty df
    intersection_df = pd.DataFrame(
        index=range(0, len(exp_df_pivot.columns) **
                    2 - len(exp_df_pivot.columns)),
        columns=["sample1", "sample2", "intersection", "fraction"])

    # populate the df
    n = 0
    for col1, col2 in itertools.combinations_with_replacement(exp_df_pivot.columns, 2):
        s1_genes = top_genes[col1]
        s2_genes = top_genes[col2]
        intersection = set(s1_genes).intersection(set(s2_genes))
        fraction = len(intersection) / float(top_n)
        intersection_df.ix[n] = [col1[1], col2[1], len(intersection), fraction]
        n += 1

        # if the samples are different, calculate the reverse intersection too
        if col1 != col2:
            intersection_df.ix[n] = [col2[1], col1[1],
                                     len(intersection), fraction]
            n += 1

    # pivot to format for heatmap.3 plotting
    intersection_df['fraction'] = intersection_df['fraction'].astype('float')
    intersection_pivot = pd.pivot_table(
        intersection_df, index="sample1", columns="sample2", values="fraction")

    for factor in set(factors_df['factor'].tolist()):
        print(factor)
        print(factors_df)
        # don't want to plot coloured by genome
        if factor == "genome":
            continue

        plotfile = "%s_%s.png" % (P.snip(outfile, ".sentinel"), factor)

        plotHeatmap = R('''
        function(int_df, fact_df){

        library(GMD)
        library(RColorBrewer)

        # subset fact_df to required factor and
        # refactor to remove unwanted levels
        fact_df = fact_df[fact_df$factor=="%(factor)s",]
        rownames(fact_df) = fact_df$sample_name
        fact_df$factor_value = factor(fact_df$factor_value)

        # set up heatmap side colours
        colours = colorRampPalette(
          brewer.pal(length(levels(fact_df$factor_value)),"Dark2"))(
            length(levels(fact_df$factor_value)))
        side_colours = colours[as.numeric((fact_df$factor_value))]
        print(side_colours)
        # plot
        png("%(plotfile)s", width=1000, heigh=1000)
        heatmap.3(as.dist(1- as.matrix(int_df)),
                  Rowv=FALSE, Colv=FALSE,
                  ColIndividualColors = side_colours,
                  RowIndividualColors = side_colours,
                  breaks=100, main="%(factor)s")
        dev.off()
        }
        ''' % locals())

        plotHeatmap(pandas2ri.py2ri(intersection_pivot),
                    pandas2ri.py2ri(factors_df))

    P.touch(outfile)


###################################################################
# Plot expression distribution
###################################################################

@mkdir("sailfish.dir/plots.dir")
@follows(loadMetaInformation,
         loadSailfishResults)
@originate("sailfish.dir/plots.dir/expression_distribution.sentinel")
def plotExpression(outfile):
    "Plot the per sample expression distibutions coloured by factor levels"

    # Note: this was being done within the pipeline but the size of
    # the dataframe seemed to be causing errors:
    # "Data values must be of type string or None."
    # See RnaseqqcReport.ExpressionDistribution tracker

    dbh = connect()

    statement = """
    SELECT sample_id, transcript_id, TPM
    FROM sailfish_transcripts"""

    df = pd.read_sql(statement, dbh)

    df['logTPM'] = df['TPM'].apply(lambda x: np.log2(x + 0.1))

    factors = dbh.execute("SELECT DISTINCT factor FROM factors")
    factors = [x[0] for x in factors if x[0] != "genome"]

    for factor in factors:

        plotfile = P.snip(outfile, ".sentinel") + "_%s.png" % factor

        factor_statement = '''
        select *
        FROM factors
        JOIN samples
        ON factors.sample_id = samples.id
        WHERE factor = "%(factor)s"''' % locals()

        factor_df = pd.read_sql(factor_statement, dbh)

        full_df = pd.merge(df, factor_df, left_on="sample_id",
                           right_on="sample_id")

        plotDistribution = R('''
        function(df){

        library(ggplot2)
        p = ggplot(df, aes(x=logTPM, group=sample_name,
                           colour=as.factor(factor_value))) +
        geom_density() +
        xlab("Log2(TPM)") + ylab("Density") +
        scale_colour_discrete(name="Factor Level") +
        theme_bw() +
        ggtitle("%(factor)s")

        ggsave("%(plotfile)s")
        }
        ''' % locals())

        plotDistribution(pandas2ri.py2ri(full_df))

    P.touch(outfile)

###################################################################
# Main pipeline tasks
###################################################################


@follows(loadContextStats,
         loadBAMStats,
         loadTranscriptProfiles,
         loadSailfishResults,
         loadMetaInformation,
         loadBias,
         loadPicardRnaSeqMetrics,
         loadAltContextStats,
         plotSailfishSaturation,
         plotTopGenesHeatmap,
         plotExpression)
def full():
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
