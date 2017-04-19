"""=====================
Read mapping pipeline
=====================

:Release: $Id$
:Date: |today|
:Tags: Python

The read mapping pipeline imports unmapped reads from one or more
NGS experiments and maps reads against a reference genome.

This pipeline works on a single genome.

Overview
========

The pipeline implements various mappers and QC plots. It can be used for

* Mapping against a genome
* Mapping RNASEQ data against a genome
* Mapping against a transcriptome

Principal targets
-----------------

mapping
    perform all mappings

qc
    perform all QC steps

full
    compute all mappings and QC

Optional targets
----------------

merge
    merge mapped :term:`bam` formatted files, for example if reads
    from different lanes were mapped separately. After merging, the
    ``qc`` target can be run again to get qc stats for the merged
    :term:`bam` formatted files.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   filename.<suffix>

The ``suffix`` determines the file type. The following suffixes/file
types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be
   sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
-------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+---------+------------+------------------------------------------------+
|*Program*|*Version*   |*Purpose*                                       |
+---------+------------+------------------------------------------------+
|bowtie_  |>=0.12.7    |read mapping                                    |
+---------+------------+------------------------------------------------+
|tophat_  |>=1.4.0     |read mapping                                    |
+---------+------------+------------------------------------------------+
|gsnap_   |>=2012.07.20|read mapping                                    |
+---------+------------+------------------------------------------------+
|samtools |>=0.1.16    |bam/sam files                                   |
+---------+------------+------------------------------------------------+
|bedtools |            |working with intervals                          |
+---------+------------+------------------------------------------------+
|sra-tools|            |extracting reads from .sra files                |
+---------+------------+------------------------------------------------+
|picard   |>=1.42      |bam/sam files. The .jar files need to be in your|
|         |            | CLASSPATH environment variable.                |
+---------+------------+------------------------------------------------+
|star_    |>=2.2.0c    |read mapping                                    |
+---------+------------+------------------------------------------------+
|bamstats_|>=1.22      |from CGR, Liverpool                             |
+---------+------------+------------------------------------------------+
|butter   |>=0.3.2     |read mapping                                    |
+---------+------------+------------------------------------------------+
|hisat    |>0.1.5      |read mapping                                    |
+---------+------------+------------------------------------------------+
|shortstack|>3.4       |read mapping                                    |
+---------+------------+------------------------------------------------+


Merging bam files
-----------------

The pipeline has the ability to merge data post-mapping. This is
useful if data have been split over several lanes and have been
provide as separate fastq files.

To enable merging, set regular expression for the input and output in
the [merge] section of the configuration file.

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz.  To run
the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz
   tar -xvzf pipeline_mapping.tgz
   cd pipeline_mapping
   python <srcdir>/pipeline_mapping.py make full

.. note::
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   tophat
      tophat_ - a read mapper to detect splice-junctions

   hisat
     hisat_ - a read mapper for RNASEQ data (basis for tophat3)

   bowtie
      bowtie_ - a read mapper

   star
      star_ - a read mapper for RNASEQ data

   bismark
      bismark_ - a read mapper for RRBS data

   butter
      butter_ - a read mapper for small RNA data (bowtie wrapper)

   shortstack - a read mapper for small RNA data (bowtie wrapper)
                that is an improvement on butter

.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _gsnap: http://research-pub.gene.com/gmap/
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011
.. _star: http://code.google.com/p/rna-star/
.. _bismark: http://www.bioinformatics.babraham.ac.uk/projects/bismark/
.. _butter: https://github.com/MikeAxtell/butter
.. _hisat: http://ccb.jhu.edu/software/hisat/manual.shtml
.. _shortstack: https://github.com/MikeAxtell/ShortStack


Code
====

"""

# load modules
from ruffus import *

import sys
import os
import re
import glob
import sqlite3
import collections

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.PipelineWindows as PipelineWindows

# Pipeline configuration
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

# Add parameters from the annotation pipeline, but
# only the interface
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))


PipelineGeneset.PARAMS = PARAMS
PipelineMappingQC.PARAMS = PARAMS

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.


# Global flags
MAPPERS = P.asList(PARAMS["mappers"])
SPLICED_MAPPING = ("tophat" in MAPPERS or
                   "gsnap" in MAPPERS or
                   "star" in MAPPERS or
                   "tophat2" in MAPPERS or
                   "transcriptome" in MAPPERS or
                   "hisat" in MAPPERS)


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


@active_if(SPLICED_MAPPING)
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

    if len(PARAMS['geneset_remove_contigs']) == 0:
        geneset_remove_contigs = None
    else:
        geneset_remove_contigs = PARAMS['geneset_remove_contigs']
    gene_ids = PipelineMapping.mergeAndFilterGTF(
        infile,
        tmp_mergedfiltered,
        "%s.removed.gz" % outfile,
        genome=os.path.join(PARAMS["genome_dir"], PARAMS["genome"]),
        max_intron_size=PARAMS["max_intron_size"],
        remove_contigs=geneset_remove_contigs,
        rna_file=rna_file)

    # Add tss_id and p_id
    PipelineMapping.resetGTFAttributes(
        infile=tmp_mergedfiltered,
        genome=os.path.join(PARAMS["genome_dir"], PARAMS["genome"]),
        gene_ids=gene_ids,
        outfile=outfile)

    os.unlink(tmp_mergedfiltered)


@active_if(SPLICED_MAPPING)
@originate("protein_coding_gene_ids.tsv")
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


@active_if(SPLICED_MAPPING)
@transform(buildReferenceGeneSet,
           suffix("reference.gtf.gz"),
           "refflat.txt")
def buildRefFlat(infile, outfile):
    '''build flat geneset for Picard RnaSeqMetrics.
    '''

    tmpflat = P.getTempFilename(".")

    statement = '''
    gtfToGenePred -genePredExt -geneNameAsName2 %(infile)s %(tmpflat)s;
    paste <(cut -f 12 %(tmpflat)s) <(cut -f 1-10 %(tmpflat)s)
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmpflat)


@active_if(SPLICED_MAPPING)
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

#########################################################################
#########################################################################
#########################################################################


@active_if(SPLICED_MAPPING)
@follows(mkdir("geneset.dir"))
@transform(PARAMS["annotations_interface_geneset_flat_gtf"],
           regex(".*"),
           add_inputs(identifyProteinCodingGenes),
           "geneset.dir/introns.gtf.gz")
def buildIntronGeneModels(infiles, outfile):
    '''build protein-coding intron-transcipts

    Retain the protein coding genes from the input gene set and
    convert the exonic sequences to intronic sequences. 10 bp is
    truncated on either end of an intron and need to have a minimum
    length of 100. Introns from nested genes might overlap, but all
    exons are removed.

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`gtf` format
    infiles[1] : str
       Input filename in :term:`tsv` format

    outfile: str
       Output filename in :term:`gtf` format

    annotations_interface_geneset_exons_gtf: str, PARAMS
       Filename for :term:`gtf` format file containing gene set exons

    '''

    filename_exons = PARAMS["annotations_interface_geneset_exons_gtf"]

    infile, genes_tsv = infiles

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter
    --map-tsv-file=%(genes_tsv)s
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort
    --sort-order=gene
    | cgat gtf2gtf
    --method=exons2introns
    --intron-min-length=100
    --intron-border=10
    --log=%(outfile)s.log
    | cgat gff2gff
    --method=crop
    --crop-gff-file=%(filename_exons)s
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=set-transcript-to-gene
    --log=%(outfile)s.log
    | awk -v OFS="\\t" -v FS="\\t" '{$3="exon"; print}'
    | gzip
    > %(outfile)s
    '''
    P.run()


@P.add_doc(PipelineGeneset.loadTranscript2Gene)
@active_if(SPLICED_MAPPING)
@transform(buildCodingGeneSet,
           suffix(".gtf.gz"),
           "_transcript2gene.load")
def loadGeneInformation(infile, outfile):
    PipelineGeneset.loadTranscript2Gene(infile, outfile)


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


@active_if(SPLICED_MAPPING)
@transform(buildCodingGeneSet, suffix(".gtf.gz"), ".fa")
def buildReferenceTranscriptome(infile, outfile):
    '''build reference transcriptome.

    Extract the sequence for each transcript in a reference geneset
    :term:`gtf` file from an indexed genome :term:`fasta` file and
    output to a :term:`fasta` file. Transcript sequences include both
    UTR and CDS.

    Additionally build :term:`bowtie` indices for tophat/tophat2 as required.

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile: str
       Output filename in :term:`fasta` format
    genome_dir : str
       :term:`PARAMS`. Directory of :term:fasta formatted files
    genome : str
       :term:`PARAMS`. Genome name (e.g hg38)
    '''

    gtf_file = P.snip(infile, ".gz")

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat %(infile)s
    | awk '$3 == "exon"' > %(gtf_file)s;
    gtf_to_fasta %(gtf_file)s %(genome_file)s %(outfile)s;
    checkpoint;
    samtools faidx %(outfile)s
    '''
    P.run()

    dest = P.snip(os.path.abspath(gtf_file), ".gtf") + ".gff"
    if not os.path.exists(dest):
        os.symlink(os.path.abspath(gtf_file), dest)

    prefix = P.snip(outfile, ".fa")

    if 'tophat' in MAPPERS or "transcriptome" in MAPPERS:
        # build raw index
        statement = '''
        bowtie-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
        '''
        P.run()

        # build color space index - disabled
        # statement = '''
        # bowtie-build -C -f %(outfile)s %(prefix)s_cs
        # >> %(outfile)s.log 2>&1
        # '''
        # P.run()

    if 'tophat2' in MAPPERS:
        statement = '''
        bowtie2-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
        '''
        P.run()

#########################################################################
#########################################################################
#########################################################################


@active_if(SPLICED_MAPPING)
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


@active_if(SPLICED_MAPPING)
@follows(mkdir("gsnap.dir"))
@merge(PARAMS["annotations_interface_geneset_exons_gtf"],
       "gsnap.dir/splicesites.iit")
def buildGSNAPSpliceSites(infile, outfile):
    '''build file with known splice sites for GSNAP from all exons

    Identify the splice from a gene set :term:`gtf` file using the
    GSNAP subprogram gts_splicesites.

    Parameters
    ----------
    infile : str
       Input filename in :term:`gtf` format
    outfile: str
       Output filename

    '''

    outfile = P.snip(outfile, ".iit")
    statement = '''zcat %(infile)s
    | gtf_splicesites | iit_store -o %(outfile)s
    > %(outfile)s.log
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
# Read mapping
#########################################################################

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    "*.remote",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|fastq.gz|fa.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz|remote)")

###################################################################
###################################################################
###################################################################
# load number of reads
###################################################################


@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\1.nreads")
def countReads(infile, outfile):
    '''Count number of reads in input files.'''
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Map reads with tophat
#########################################################################


@active_if(SPLICED_MAPPING)
@follows(mkdir("tophat.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildJunctions, buildReferenceTranscriptome),
           r"tophat.dir/\1.tophat.bam")
def mapReadsWithTophat(infiles, outfile):
    """
    Map reads using Tophat (spliced reads).

    Parameters
    ----------

    infiles: list
        contains 3 filenames -
    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        :term:`fasta` filename, suffix .fa
        reference transcriptome

    infiles[2]: str
        filename with suffix .junctions containing a list of known
        splice junctions.

    tophat_threads: int
        :term:`PARAMS`
        number of threads with which to run tophat

    tophat_options: str
        :term:`PARAMS`
        string containing options to pass to tophat

    tophat_memory: str
        :term:`PARAMS`
        memory required for tophat job

    tophat_executable: str
        :term:`PARAMS`
        path to tophat executable

    strandness
        :term:`PARAMS`
        FR, RF, F or R or empty see
        http://www.ccb.jhu.edu/software/hisat/manual.shtml#options
        will be converted to tophat specific option

    tophat_include_reference_transcriptome: bool
        :term:`PARAMS`
        if set, map to reference transcriptome

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    bowtie_index_dir: str
        :term:`PARAMS`
        path to directory containing bowtie indices

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.


    .. note::
       If tophat fails with an error such as::

          Error: segment-based junction search failed with err =-6
          what():  std::bad_alloc

       it means that it ran out of memory.

    """

    job_threads = PARAMS["tophat_threads"]

    # convert strandness to tophat-style library type
    if PARAMS["strandness"] == ("RF" or "R"):
        tophat_library_type = "fr-firststrand"
    elif PARAMS["strandness"] == ("FR" or "F"):
        tophat_library_type = "fr-secondstrand"
    else:
        tophat_library_type = "fr-unstranded"

    if "--butterfly-search" in PARAMS["tophat_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_memory = "50G"
    else:
        job_memory = PARAMS["tophat_memory"]

    m = PipelineMapping.Tophat(
        executable=P.substituteParameters(**locals())["tophat_executable"],
        strip_sequence=PARAMS["strip_sequence"],
        tool_options=PARAMS["tophat_options"])
    infile, reffile, transcriptfile = infiles
    tophat_options = PARAMS["tophat_options"] + \
        " --raw-juncs %(reffile)s " % locals()

    # Nick - added the option to map to the reference transcriptome first
    # (built within the pipeline)
    if PARAMS["tophat_include_reference_transcriptome"]:
        prefix = os.path.abspath(P.snip(transcriptfile, ".fa"))
        tophat_options = tophat_options + \
            " --transcriptome-index=%s -n 2" % prefix

    statement = m.build((infile,), outfile)
    P.run()


@active_if(SPLICED_MAPPING)
@follows(mkdir("tophat2.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildJunctions, buildReferenceTranscriptome),
           r"tophat2.dir/\1.tophat2.bam")
def mapReadsWithTophat2(infiles, outfile):
    '''
     Map reads using Tophat2 (spliced reads).

    Parameters
    ----------

    infiles: list
        contains 3 filenames -
    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        :term:`fasta` filename, suffix .fa
        reference transcriptome

    infiles[2]: str
        filename with suffix .junctions containing a list of known
        splice junctions.

    tophat2_threads: int
        :term:`PARAMS`
        number of threads with which to run tophat2

    tophat2_options: str
        :term:`PARAMS`
        string containing options to pass to tophat2

    tophat2_memory: str
        :term:`PARAMS`
        memory required for tophat2 job

    tophat2_executable: str
        :term:`PARAMS`
        path to tophat2 executable

    strandness
        :term:`PARAMS`
        FR, RF, F or R or empty see
        http://www.ccb.jhu.edu/software/hisat/manual.shtml#options
        will be converted to tophat specific option

    tophat2_include_reference_transcriptome: bool
        :term:`PARAMS`
        if set, map to reference transcriptome

    tophat2_mate_inner_dist: int
        :term:`PARAMS`
        insert length (2 * read length)

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    bowtie_index_dir: str
        :term:`PARAMS`
        path to directory containing bowtie indices

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.


    .. note::
       If tophat fails with an error such as::

          Error: segment-based junction search failed with err =-6
          what():  std::bad_alloc

       it means that it ran out of memory.

    '''
    job_threads = PARAMS["tophat2_threads"]

    # convert strandness to tophat-style library type
    if PARAMS["strandness"] == ("RF" or "R"):
        tophat2_library_type = "fr-firststrand"
    elif PARAMS["strandness"] == ("FR" or "F"):
        tophat2_library_type = "fr-secondstrand"
    else:
        tophat2_library_type = "fr-unstranded"

    if "--butterfly-search" in PARAMS["tophat2_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_memory = "50G"
    else:
        job_memory = PARAMS["tophat2_memory"]

    m = PipelineMapping.Tophat2(
        executable=P.substituteParameters(**locals())["tophat2_executable"],
        strip_sequence=PARAMS["strip_sequence"],
        tool_options=PARAMS["tophat2_options"])

    infile, reffile, transcriptfile = infiles
    tophat2_options = PARAMS["tophat2_options"] + \
        " --raw-juncs %(reffile)s" % locals()

    # Nick - added the option to map to the reference transcriptome first
    # (built within the pipeline)
    if PARAMS["tophat2_include_reference_transcriptome"]:
        prefix = os.path.abspath(P.snip(transcriptfile, ".fa"))
        tophat2_options = tophat2_options + \
            " --transcriptome-index=%s -n 2" % prefix

    statement = m.build((infile,), outfile)
    P.run()

############################################################
############################################################
############################################################


@active_if(SPLICED_MAPPING)
@follows(mkdir("hisat.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
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

    strandness: str
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
        executable=P.substituteParameters(**locals())["hisat_executable"],
        strip_sequence=PARAMS["strip_sequence"],
        stranded=PARAMS["strandness"])

    infile, junctions = infiles

    statement = m.build((infile,), outfile)

    P.run()

############################################################
############################################################
############################################################


@active_if(SPLICED_MAPPING)
@merge(mapReadsWithTophat, "tophat_stats.tsv")
def buildTophatStats(infiles, outfile):
    '''
    Build stats about tophat runs.

    Uses the log files from tophat mapping runs to build a table showing
    counts for various statistics for each input file.
    These statistics are: reads in, reads removed, reads out, junctions loaded,
    junctions found, possible splices.

    Parameters
    ----------
    infiles: list
        list of filenames of :term:`bam` files containing mapped reads

    paired_end: bool
        :term:`PARAMS` if true, reads are paired end.

    outfile: str
        :term:`tsv` file to write the stats about the run

    '''
    def _select(lines, pattern):
        '''
        Looks for a pattern in each line of the bam file

        Parameters
        ----------
        lines: list
            readlines object from log file from the tophat run
        pattern: str
            regex specifying the pattern to search for.
        '''
        x = re.compile(pattern)
        for line in lines:
            r = x.search(line)
            if r:
                g = r.groups()
                if len(g) > 1:
                    return g
                else:
                    return g[0]

        raise ValueError("pattern '%s' not found %s" % (pattern, lines))

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join(("track",
                          "reads_in",
                          "reads_removed",
                          "reads_out",
                          "junctions_loaded",
                          "junctions_found",
                          "possible_splices")) + "\n")

    for infile in infiles:

        track = P.snip(infile, ".bam")
        indir = infile + ".logs"

        fn = os.path.join(indir, "prep_reads.log")
        lines = open(fn).readlines()
        reads_removed, reads_in = list(map(
            int, _select(lines, "(\d+) out of (\d+) reads have been filtered out")))
        reads_out = reads_in - reads_removed
        prep_reads_version = _select(lines, "prep_reads (.*)$")

        fn = os.path.join(indir, "reports.log")
        lines = open(fn).readlines()
        tophat_reports_version = _select(lines, "tophat_reports (.*)$")
        junctions_loaded = int(_select(lines, "Loaded (\d+) junctions"))
        junctions_found = int(
            _select(lines, "Found (\d+) junctions from happy spliced reads"))

        fn = os.path.join(indir, "segment_juncs.log")

        if os.path.exists(fn):
            lines = open(fn).readlines()
            if len(lines) > 0:
                segment_juncs_version = _select(lines, "segment_juncs (.*)$")
                possible_splices = int(
                    _select(lines, "Reported (\d+) total potential splices"))
            else:
                segment_juncs_version = "na"
                possible_splices = ""
        else:
            segment_juncs_version = "na"
            possible_splices = ""

        # fix for paired end reads - tophat reports pairs, not reads
        if PARAMS["paired_end"]:
            reads_in *= 2
            reads_out *= 2
            reads_removed *= 2

        outf.write("\t".join(map(str, (
            track,
            reads_in, reads_removed, reads_out,
            junctions_loaded, junctions_found, possible_splices))) + "\n")

    outf.close()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@transform(buildTophatStats, suffix(".tsv"), ".load")
def loadTophatStats(infile, outfile):
    '''
    Loads statistics about a tophat run from a tsv file to a database table -
    tophat_stats.
    The columns are track, reads in, reads removed, reads out, junctions loaded
    , junctions found, possible splices.

    Parameters
    ----------
    infile: term:`tsv` file containing a table of tophat statistics.
    outfile: .load file
    '''
    P.load(infile, outfile)


@active_if(SPLICED_MAPPING)
@follows(mkdir("gsnap.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildGSNAPSpliceSites),
           r"gsnap.dir/\1.gsnap.bam")
def mapReadsWithGSNAP(infiles, outfile):
    '''
    Maps reads using GSNAP (mRNA and EST sequences).

    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        filename of type iit containing all known splice sites

    gsnap_memory: str
        :term:`PARAMS`
        memory required for gsnap job

    gsnap_node_threads: int
        :term:`PARAMS`
        number of threads to use on the node

    gsnap_worker_threads: int
        :term:`PARAMS`
        --nthreads option for GSNAP described in
        http://research-pub.gene.com/gmap/src/README
        this number of threads plus 2 will be used

    gsnap_options: str
        :term:`PARAMS`
        string containing command line options for GSNAP,
        details at http://research-pub.gene.com/gmap/src/README

    gsnap_executable: str
        :term:`PARAMS`
        path to gsnap executable

    gsnap_mapping_genome: str
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    infile, infile_splices = infiles

    job_memory = PARAMS["gsnap_memory"]
    job_threads = PARAMS["gsnap_node_threads"]
    gsnap_mapping_genome = PARAMS["gsnap_genome"] or PARAMS["genome"]

    m = PipelineMapping.GSNAP(
        executable=P.substituteParameters(**locals())["gsnap_executable"],
        strip_sequence=PARAMS["strip_sequence"])

    if PARAMS["gsnap_include_known_splice_sites"]:
        gsnap_options = PARAMS["gsnap_options"] + \
            " --use-splicing=%(infile_splices)s " % locals()

    statement = m.build((infile,), outfile)
    P.run()


@active_if(SPLICED_MAPPING)
@follows(mkdir("star.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"star.dir/\1.star.bam")
def mapReadsWithSTAR(infile, outfile):
    '''
    Maps reads using STAR (spliced reads).

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    star_memory: str
        :term:`PARAMS`
        memory required for STAR job

    star_threads: int
        :term:`PARAMS`
        number of threads with which to run STAR

    star_genome: str
        :term:`PARAMS`
        path to genome if using a splice junction database sjdb

    genome: str
        :term:`PARAMS`
        path to genome if not using a splice junction database

    star_executable: str
        :term:`PARAMS`
        path to star executable

    star_index_dir: str
        :term:`PARAMS`
        path to directory containing star indices.

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.

    '''

    job_threads = PARAMS["star_threads"]
    job_memory = PARAMS["star_memory"]

    star_mapping_genome = PARAMS["star_genome"] or PARAMS["genome"]

    m = PipelineMapping.STAR(
        executable=P.substituteParameters(**locals())["star_executable"],
        strip_sequence=PARAMS["strip_sequence"])

    statement = m.build((infile,), outfile)
    P.run()


@active_if(SPLICED_MAPPING)
@merge(mapReadsWithSTAR, "star_stats.tsv")
def buildSTARStats(infiles, outfile):
    '''Compile statistics from STAR run

    Concatenates log files from STAR runs and reformats them as a tab
    delimited table.

    Parameters
    ----------
    infile: list
        :term:`bam` files generated with STAR.
    outfile: str
        :term: `tsv` file containing statistics about STAR run
    '''

    data = collections.defaultdict(list)
    for infile in infiles:
        fn = infile + ".final.log"
        if not os.path.exists(fn):
            raise ValueError("incomplete run: %s" % infile)

        for line in IOTools.openFile(fn):
            if "|" not in line:
                continue
            header, value = line.split("|")
            header = re.sub("%", "percent", header)
            data[header.strip()].append(value.strip())

    keys = list(data.keys())
    outf = IOTools.openFile(outfile, "w")
    outf.write("track\t%s\n" % "\t".join(keys))
    for x, infile in enumerate(infiles):
        track = P.snip(os.path.basename(infile), ".bam")
        outf.write("%s\t%s\n" %
                   (track, "\t".join([data[key][x] for key in keys])))
    outf.close()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@transform(buildSTARStats, suffix(".tsv"), ".load")
def loadSTARStats(infile, outfile):
    '''
    Loads statistics about a star run from a tsv file to a database table -
    star_stats.

    Parameters
    ----------
    infile: term:`tsv` file containing a table of tophat statistics.
    outfile: .load file logging database loading
    '''
    P.load(infile, outfile)


@active_if(SPLICED_MAPPING)
@follows(mkdir("transcriptome.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(buildReferenceTranscriptome),
           r"transcriptome.dir/\1.trans.bam")
def mapReadsWithBowtieAgainstTranscriptome(infiles, outfile):
    '''
    Map reads using bowtie against transcriptome data.

    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        :term:`fasta` file containing reference genome

    bowtie_threads: int
        :term:`PARAMS`
        number of threads with which to run bowtie

    bowtie_memory: str
        :term:`PARAMS`
        memory required for bowtie job

    bowtie_executable: str
        :term:`PARAMS`
        path to bowtie executable

    bowtie_index_dir: str
        :term:`PARAMS`
        path to directory containing bowtie indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.

    '''

    # Mapping will permit up to one mismatches. This is sufficient
    # as the downstream filter in rnaseq_bams2bam requires the
    # number of mismatches less than the genomic number of mismatches.
    # Change this, if the number of permitted mismatches for the genome
    # increases.

    # Output all valid matches in the best stratum. This will
    # inflate the file sizes due to matches to alternative transcripts
    # but otherwise matches to paralogs will be missed (and such
    # reads would be filtered out).
    job_threads = PARAMS["bowtie_threads"]
    m = PipelineMapping.BowtieTranscripts(
        executable=P.substituteParameters(**locals())["bowtie_executable"],
        strip_sequence=PARAMS["strip_sequence"])
    infile, reffile = infiles
    prefix = P.snip(reffile, ".fa")
    # IMS: moved reporting options to ini
    # bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_transcriptome_options"]
    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("bowtie.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(
               os.path.join(PARAMS["bowtie_index_dir"],
                            PARAMS["genome"] + ".fa")),
           r"bowtie.dir/\1.bowtie.bam")
def mapReadsWithBowtie(infiles, outfile):
    '''
    Map reads with bowtie (short reads).
    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        :term:`fasta` file containing reference genome

    bowtie_threads: int
        :term:`PARAMS`
        number of threads with which to run bowtie

    bowtie_memory: str
        :term:`PARAMS`
        memory required for bowtie job

    bowtie_executable: str
        :term:`PARAMS`
        path to bowtie executable

    bowtie_options: str
        :term:`PARAMS`
        string containing command line options for bowtie - refer
        to http://bowtie-bio.sourceforge.net/index.shtml

    bowtie_index_dir: str
        :term:`PARAMS`
        path to directory containing bowtie indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    job_threads = PARAMS["bowtie_threads"]
    job_memory = PARAMS["bowtie_memory"]

    m = PipelineMapping.Bowtie(
        executable=P.substituteParameters(**locals())["bowtie_executable"],
        tool_options=P.substituteParameters(**locals())["bowtie_options"],
        strip_sequence=PARAMS["strip_sequence"])
    infile, reffile = infiles
    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("bowtie2.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           add_inputs(
               os.path.join(PARAMS["bowtie_index_dir"],
                            PARAMS["genome"] + ".fa")),
           r"bowtie2.dir/\1.bowtie2.bam")
def mapReadsWithBowtie2(infiles, outfile):
    '''
    Map reads with bowtie2.
    Parameters
    ----------
    infiles: list
        contains two filenames -

    infiles[0]: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    infiles[1]: str
        :term:`fasta` file containing reference genome

    bowtie2_threads: int
        :term:`PARAMS`
        number of threads with which to run bowtie2

    bowtie2_memory: str
        :term:`PARAMS`
        memory required for bowtie2 job

    bowtie2_executable: str
        :term:`PARAMS`
        path to bowtie2 executable

    bowtie2_options: str
        :term:`PARAMS`
        string containing command line options for bowtie2 -
        refer to http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    job_threads = PARAMS["bowtie2_threads"]
    job_memory = PARAMS["bowtie2_memory"]

    m = PipelineMapping.Bowtie2(
        executable=P.substituteParameters(**locals())["bowtie2_executable"],
        tool_options=P.substituteParameters(**locals())["bowtie2_options"],
        strip_sequence=PARAMS["strip_sequence"])
    infile, reffile = infiles
    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("bwa.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"bwa.dir/\1.bwa.bam")
def mapReadsWithBWA(infile, outfile):
    '''
    Map reads with bwa

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    bwa_threads: int
        :term:`PARAMS`
        number of threads with which to run BWA

    bwa_memory: str
        :term:`PARAMS`
        memory required for BWA job

    bwa_algorithm: str
        :term:`PARAMS`
        two options - 'aln' or 'mem' - refer to
        http://bio-bwa.sourceforge.net/bwa.shtml

    bwa_set_nh: str
        :term:`PARAMS`
        sets the NH tag    which specifies multiple hits in sam format
        otherwise not set by bwa

    bwa_index_dir: str
        :term:`PARAMS`
        path to directory containing bwa indices

    bwa_aln_options: str
        :term:`PARAMS`
        string containing parameters for bwa if run using the
        'aln' algorithm - refer to
        http://bio-bwa.sourceforge.net/bwa.shtml

    bwa_samse_options: str
        :term:`PARAMS`
        string containing single end read options for bwa - refer
        to http://bio-bwa.sourceforge.net/bwa.shtml

    bwa_mem_options: str
        :term:`PARAMS`
        string containing parameters for bwa if run using the
        'mem' algorithm    - refer to
        http://bio-bwa.sourceforge.net/bwa.shtml

    remove_non_unique
        :term:`PARAMS`
        If true, a filtering step is included in postprocess, which removes
        reads that have more than 1 best hit

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    job_threads = PARAMS["bwa_threads"]
    job_memory = PARAMS["bwa_memory"]

    if PARAMS["bwa_algorithm"] == "aln":
        m = PipelineMapping.BWA(
            remove_non_unique=PARAMS["remove_non_unique"],
            strip_sequence=PARAMS["strip_sequence"],
            set_nh=PARAMS["bwa_set_nh"])
    elif PARAMS["bwa_algorithm"] == "mem":
        m = PipelineMapping.BWAMEM(
            remove_non_unique=PARAMS["remove_non_unique"],
            strip_sequence=PARAMS["strip_sequence"],
            set_nh=PARAMS["bwa_set_nh"])
    else:
        raise ValueError("bwa algorithm '%s' not known" % algorithm)

    statement = m.build((infile,), outfile)
    P.run()


@follows(mkdir("stampy.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"stampy.dir/\1.stampy.bam")
def mapReadsWithStampy(infile, outfile):
    '''
    Map reads with stampy

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    genome
        :term:`PARAMS`
        path to reference genome

    stampy_threads: int
        :term:`PARAMS`
        number of threads with which to run Stampy

    stampy_memory: str
        :term:`PARAMS`
        memory required for stampy job

    stampy_index_dir: str
        :term:`PARAMS`
        path to directory containing stampy indices

    bwa_index_dir: str
        :term:`PARAMS`
        path to directory containing bwa indices

    strip_sequence: bool
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    job_threads = PARAMS["stampy_threads"]
    job_memory = PARAMS["stampy_memory"]

    m = PipelineMapping.Stampy(strip_sequence=PARAMS["strip_sequence"])
    statement = m.build((infile,), outfile)
    P.run()

###################################################################
###################################################################
###################################################################
# Map reads with butter
###################################################################


@follows(mkdir("butter.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"butter.dir/\1.butter.bam")
def mapReadsWithButter(infile, outfile):
    '''
    Map reads with butter

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    butter_threads: int
        :term:`PARAMS`
        number of threads with which to run butter

    butter_memory: str
        :term:`PARAMS`
        memory required for butter job

    butter_options: str
        :term:`PARAMS`
        string containing command line options to pass to Butter -
        refer to https://github.com/MikeAxtell/butter

    butter_index_dir: str
        :term:`PARAMS`
        path to directory containing butter indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    butter_set_nh: str
        :term:`PARAMS`
        sets the NH tag    which specifies multiple hits in sam format
        otherwise not set by butter

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    # easier to check whether infiles are paired reads here
    if infile.endswith(".sra"):
        outdir = P.getTempDir()
        f = Sra.sneak(infile, outdir)
        shutil.rmtree(outdir)
        assert len(f) == 1, NotImplementedError('''The sra archive contains
        paired end data,Butter does not support paired end reads''')

    elif infile.endswith(".csfasta.F3.gz") or infile.endswith(".fastq.1.gz"):
        raise NotImplementedError('''infiles are paired end: %(infile)s,
        Butter does not support paired end reads''' % locals())

    job_threads = PARAMS["butter_threads"]
    job_memory = PARAMS["butter_memory"]

    m = PipelineMapping.Butter(
        strip_sequence=PARAMS["strip_sequence"],
        set_nh=PARAMS["butter_set_nh"])
    statement = m.build((infile,), outfile)

    P.run()

###################################################################
###################################################################
###################################################################
# Map reads with shortstack
###################################################################


@follows(mkdir("shortstack.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"shortstack.dir/\1.shortstack.bam")
def mapReadsWithShortstack(infile, outfile):
    '''
    Map reads with shortstack

    Parameters
    ----------
    infile: str
        filename of reads file
        can be :term:`fastq`, :term:`sra`, csfasta

    shortstack_threads: int
        :term:`PARAMS`
        number of threads with which to run shortstack

    shortstack_memory: str
        :term:`PARAMS`
        memory required for shortstack job

    shortstack_options: str
        :term:`PARAMS`
        string containing command line options to pass to Shortstack -
        refer to https://github.com/MikeAxtell/ShortStack

    shortstack_index_dir: str
        :term:`PARAMS`
        path to directory containing shortstack indices

    strip_sequence: bool
        :term:`PARAMS`
        if set, strip read sequence and quality information

    outfile: str
        :term:`bam` filename to write the mapped reads in bam format.
    '''

    # easier to check whether infiles are paired reads here
    if infile.endswith(".sra"):
        outdir = P.getTempDir()
        f = Sra.sneak(infile, outdir)
        shutil.rmtree(outdir)
        assert len(f) == 1, NotImplementedError('''The sra archive contains
        paired end data,shortstack does not support paired end reads''')

    elif infile.endswith(".csfasta.F3.gz") or infile.endswith(".fastq.1.gz"):
        raise NotImplementedError('''infiles are paired end: %(infile)s,
        shortstack does not support paired end reads''' % locals())

    job_threads = PARAMS["shortstack_threads"]
    job_memory = PARAMS["shortstack_memory"]

    m = PipelineMapping.Shortstack(
        strip_sequence=PARAMS["strip_sequence"])
    statement = m.build((infile,), outfile)

    P.run()

###################################################################
###################################################################
###################################################################
# Create map reads tasks
###################################################################


MAPPINGTARGETS = []
mapToMappingTargets = {'tophat': (mapReadsWithTophat, loadTophatStats),
                       'tophat2': (mapReadsWithTophat2,),
                       'bowtie': (mapReadsWithBowtie,),
                       'bowtie2': (mapReadsWithBowtie2,),
                       'bwa': (mapReadsWithBWA,),
                       'stampy': (mapReadsWithStampy,),
                       'transcriptome':
                       (mapReadsWithBowtieAgainstTranscriptome,),
                       'gsnap': (mapReadsWithGSNAP,),
                       'star': (mapReadsWithSTAR, loadSTARStats),
                       'butter': (mapReadsWithButter,),
                       'shortstack': (mapReadsWithShortstack,),
                       'hisat': (mapReadsWithHisat,)
                       }

for x in P.asList(PARAMS["mappers"]):
    MAPPINGTARGETS.extend(mapToMappingTargets[x])


@follows(*MAPPINGTARGETS)
def mapping():
    ''' dummy task to define upstream mapping tasks'''
    pass


if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    if "merge_pattern_output" not in PARAMS or \
       not PARAMS["merge_pattern_output"]:
        raise ValueError(
            "no output pattern 'merge_pattern_output' specified")

    @collate(MAPPINGTARGETS,
             regex("%s\.([^.]+).bam" % PARAMS["merge_pattern_input"].strip()),
             # the last expression counts number of groups in pattern_input
             r"%s.\%i.bam" % (PARAMS["merge_pattern_output"].strip(),
                              PARAMS["merge_pattern_input"].count("(") + 1),
             )
    def mergeBAMFiles(infiles, outfile):
        '''merge BAM files from the same experiment using user-defined regex

        For the mapping stages it is beneficial to perform mapping
        seperately for each sequence read infile(s) per sample so that
        the consistency can be checked. However, for downstream tasks,
        the merged :term:`bam` alignment files are required.

        Parameters
        ----------
        infiles : list
           list of :term:`bam` format alignment files
        outfile : str
           Output filename in :term:`bam` format
        '''

        if len(infiles) == 1:
            E.info(
                "%(outfile)s: only one file for merging - creating "
                "softlink" % locals())
            P.clone(infiles[0], outfile)
            P.clone(infiles[0] + ".bai", outfile + ".bai")
            return

        infiles = " ".join(infiles)
        statement = '''
        samtools merge %(outfile)s %(infiles)s >& %(outfile)s.log;
        checkpoint;
        samtools index %(outfile)s
        '''
        P.run()

    MAPPINGTARGETS = MAPPINGTARGETS + [mergeBAMFiles]

    @collate(countReads,
             regex("%s.nreads" % PARAMS["merge_pattern_input"]),
             r"%s.nreads" % PARAMS["merge_pattern_output"],
             )
    def mergeReadCounts(infiles, outfile):
        '''merge read counts files from the same experiment using
        user-defined regex

        For the mapping stages it is beneficial to perform mapping
        seperately for each sequence read infile(s) so that the
        consistency can be checked. However, for downstream tasks, the
        merged counts per sample are required.

        Parameters
        ----------
        infiles : list of str
           list of filenames containing read counts per sequence read infile
        outfile : str
           Output filename containing total counts for a sample
        '''

        nreads = 0
        for infile in infiles:
            with IOTools.openFile(infile, "r") as inf:
                for line in inf:
                    if not line.startswith("nreads"):
                        continue
                    E.info("%s" % line[:-1])
                    nreads += int(line[:-1].split("\t")[1])

        outf = IOTools.openFile(outfile, "w")
        outf.write("nreads\t%i\n" % nreads)
        outf.close()

else:
    @follows(countReads)
    @transform(SEQUENCEFILES,
               SEQUENCEFILES_REGEX,
               r"nreads.dir/\1.nreads")
    # this decorator for the dummy mergeReadCounts is needed to prevent
    # rerunning of all downstream functions.
    def mergeReadCounts(infiles, outfiles):
        pass

###################################################################
###################################################################
###################################################################
# QC targets
###################################################################


@P.add_doc(PipelineMappingQC.buildPicardAlignmentStats)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"] + ".fa")),
           ".picard_stats")
def buildPicardStats(infiles, outfile):
    ''' build Picard alignment stats '''
    infile, reffile = infiles

    # patch for mapping against transcriptome - switch genomic reference
    # to transcriptomic sequences
    if "transcriptome.dir" in infile:
        reffile = "refcoding.fa"

    PipelineMappingQC.buildPicardAlignmentStats(infile,
                                                outfile,
                                                reffile)


@P.add_doc(PipelineMappingQC.loadPicardAlignmentStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardStats, "picard_stats.load")
def loadPicardStats(infiles, outfile):
    '''merge alignment stats into single tables.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)


@P.add_doc(PipelineMappingQC.buildPicardDuplicationStats)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           ".picard_duplication_metrics")
def buildPicardDuplicationStats(infile, outfile):
    '''Get duplicate stats from picard MarkDuplicates '''
    PipelineMappingQC.buildPicardDuplicationStats(infile, outfile)


@P.add_doc(PipelineMappingQC.loadPicardDuplicationStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildPicardDuplicationStats, ["picard_duplication_stats.load",
                                     "picard_duplication_histogram.load"])
def loadPicardDuplicationStats(infiles, outfiles):
    '''merge alignment stats into single tables.'''
    # separate load function while testing
    PipelineMappingQC.loadPicardDuplicationStats(infiles, outfiles)


@follows(countReads, mergeReadCounts)
@transform(MAPPINGTARGETS,
           regex("(.*)/(.*)\.(.*).bam"),
           add_inputs(r"nreads.dir/\2.nreads"),
           r"\1/\2.\3.readstats")
def buildBAMStats(infiles, outfile):
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

    job_memory = "32G"

    bamfile, readsfile = infiles

    nreads = PipelineMappingQC.getNumReadsFromReadsFile(readsfile)
    track = P.snip(os.path.basename(readsfile),
                   ".nreads")

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
         --num-reads=%(nreads)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(bamfile)s
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
@transform(MAPPINGTARGETS,
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


@transform(MAPPINGTARGETS,
           suffix(".bam"),
           ".idxstats")
def buildIdxstats(infile, outfile):
    '''gets idxstats for bam file so number of reads per chromosome can
    be plotted later'''

    statement = '''samtools idxstats %(infile)s > %(outfile)s''' % locals()
    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildIdxstats, "idxstats_reads_per_chromosome.load")
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

###################################################################
###################################################################
###################################################################
# QC specific to spliced mapping
###################################################################


@active_if(SPLICED_MAPPING)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(buildCodingExons),
           ".exon.validation.tsv.gz")
def buildExonValidation(infiles, outfile):
    '''Compare the alignments to the exon models to quantify exon
    overrun/underrun

    Expectation is that reads should not extend beyond known exons.

    Parameters
    ----------
    infiles : list
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`gtf` format with exon validation stats
    '''

    infile, exons = infiles
    statement = '''cat %(infile)s
    | cgat bam_vs_gtf
         --exons-file=%(exons)s
         --force-output
         --log=%(outfile)s.log
         --output-filename-pattern="%(outfile)s.%%s.gz"
    | gzip
    > %(outfile)s
    '''

    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@merge(buildExonValidation, "exon_validation.load")
def loadExonValidation(infiles, outfile):
    ''' load individual and merged exon validation stats

    For each sample, the exon validation stats are loaded into a table
    named by sample and mapper
    [sample]_[mapper]_overrun

    The merge alignment stats for all samples are merged and loaded
    into single table called exon_validation

    Parameters
    ----------
    infiles : list
       Input filenames with exon validation stats
    outfile : str
       Output filename
    '''

    suffix = ".exon.validation.tsv.gz"

    P.mergeAndLoad(infiles, outfile, suffix=suffix)
    for infile in infiles:
        track = P.snip(infile, suffix)
        o = "%s_overrun.load" % track
        P.load(infile + ".overrun.gz", o)


@active_if(SPLICED_MAPPING)
@transform(MAPPINGTARGETS,
           regex("(.+).bam"),
           add_inputs(buildCodingGeneSet),
           r"\1.transcript_counts.tsv.gz")
def buildTranscriptLevelReadCounts(infiles, outfile):
    '''count reads in gene models

    Count the reads from a :term:`bam` file which overlap the
    positions of protein coding transcripts in a :term:`gtf` format
    transcripts file.

    Parameters
    ----------
    infiles : list of str
    infiles[0] : str
       Input filename in :term:`bam` format
    infiles[1] : str
       Input filename in :term:`gtf` format

    outfile : str
       Output filename in :term:`tsv` format


    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.

    '''
    infile, geneset = infiles

    job_memory = "8G"

    statement = '''
    zcat %(geneset)s
    | cgat gtf2table
    --reporter=transcripts
    --bam-file=%(infile)s
    --counter=length
    --column-prefix="exons_"
    --counter=read-counts
    --column-prefix=""
    --counter=read-coverage
    --column-prefix=coverage_
    -v 0
    | gzip
    > %(outfile)s
    ''' % locals()

    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@transform(buildTranscriptLevelReadCounts,
           suffix(".tsv.gz"),
           ".load")
def loadTranscriptLevelReadCounts(infile, outfile):
    ''' load the transcript level read counts and index on the transcript_id

    Parameters
    ----------
    infile : str
       Input filename in :term:`tsv` format
    outfile : str
       Output filename, the table name is derived from `outfile`

    '''
    P.load(infile, outfile,
           options="--add-index=transcript_id --allow-empty-file")


@active_if(SPLICED_MAPPING)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(buildIntronGeneModels),
           ".intron_counts.tsv.gz")
def buildIntronLevelReadCounts(infiles, outfile):
    '''count reads in gene models
    Count the reads from a :term:`bam` file which overlap the
    positions of introns in a :term:`gtf` format transcripts file.
    Parameters
    ----------
    infiles : list of str
       infile :term:`str`
          Input filename in :term:`bam` format
       geneset :term:`str`
          Input filename in :term:`gtf` format
    outfile : str
       Output filename in :term:`tsv` format
    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
    '''
    infile, exons = infiles

    job_memory = "4G"

    if "transcriptome.dir" in infile:
        P.touch(outfile)
        return

    statement = '''
    zcat %(exons)s
    | cgat gtf2table
          --reporter=genes
          --bam-file=%(infile)s
          --counter=length
          --column-prefix="introns_"
          --counter=read-counts
          --column-prefix=""
          --counter=read-coverage
          --column-prefix=coverage_
    | gzip
    > %(outfile)s
    '''

    P.run()


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@active_if(SPLICED_MAPPING)
@transform(buildIntronLevelReadCounts,
           suffix(".tsv.gz"),
           ".load")
def loadIntronLevelReadCounts(infile, outfile):
    ''' load the intron level read counts and index on the transcript_id

    Parameters
    ----------
    infiles : str
       Input filename in :term:`tsv` format
    outfile : str
       Output filename, the table name is derived from `outfile`

    '''
    P.load(infile, outfile, options="--add-index=gene_id --allow-empty-file")


@merge((countReads, mergeReadCounts), "reads_summary.load")
def loadReadCounts(infiles, outfile):
    ''' load the read counts

    individual read counts are merged and loaded into a table called
    reads_summary

    Parameters
    ----------
    infile : str
       Input filename in :term:`tsv` format
    outfile : str
       Output filename, the table name is derived from `outfile`

    '''

    outf = P.getTempFile(".")
    outf.write("track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile(infile).readlines()
        nreads = int(lines[0][:-1].split("\t")[1])
        outf.write("%s\t%i\n" % (track, nreads))
    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)


@active_if(SPLICED_MAPPING)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(buildCodingExons),
           ".transcriptprofile.gz")
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
    --method=geneprofileabsolutedistancefromthreeprimeend
    --normalize-profile=all
    %(bamfile)s %(gtffile)s
    | gzip
    > %(outfile)s
    '''

    P.run()


@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(buildCodingExons),
           ".geneprofile.gz")
def buildGeneProfiles(infiles, outfile):
    '''build gene coverage profiles

    Gene coverage plots are useful to determine mapping location
    for a number of genomic techniques.

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
                --method=geneprofilewithintrons
                --bam-file=%(bamfile)s
                --gtf-file=%(gtffile)s
                --normalize-transcript=total-sum
                --normalize-profile=area
                --log=%(outfile)s.log
                --output-filename-pattern=%(outfile)s.%%s
                > %(outfile)s '''
    P.run()


@active_if(SPLICED_MAPPING)
@P.add_doc(PipelineMappingQC.buildPicardRnaSeqMetrics)
@transform(MAPPINGTARGETS,
           suffix(".bam"),
           add_inputs(buildRefFlat),
           ".picard_rna_metrics")
def buildPicardRnaSeqMetrics(infiles, outfile):
    '''Get duplicate stats from picard RNASeqMetrics '''
    # convert strandness to tophat-style library type
    if PARAMS["strandness"] == ("RF" or "R"):
        strand = "SECOND_READ_TRANSCRIPTION_STRAND"
    elif PARAMS["strandness"] == ("FR" or "F"):
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
###################################################################
###################################################################
# various export functions
###################################################################


@transform(MAPPINGTARGETS,
           regex(".bam"),
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

    if SPLICED_MAPPING:
        # use bedtools for RNASEQ data

        # scale by Million reads mapped
        reads_mapped = BamTools.getNumberOfAlignments(infile)
        scale = 1000000.0 / float(reads_mapped)
        tmpfile = P.getTempFilename()
        contig_sizes = PARAMS["annotations_interface_contigs"]
        job_memory = "3G"
        statement = '''bedtools genomecov
        -ibam %(infile)s
        -g %(contig_sizes)s
        -bg
        -split
        -scale %(scale)f
        > %(tmpfile)s;
        checkpoint;
        bedGraphToBigWig %(tmpfile)s %(contig_sizes)s %(outfile)s;
        checkpoint;
        rm -f %(tmpfile)s
        '''
    else:
        # wigToBigWig observed to use 16G
        job_memory = "16G"
        statement = '''cgat bam2wiggle
        --output-format=bigwig
        %(bigwig_options)s
        %(infile)s
        %(outfile)s
        > %(outfile)s.log'''
    P.run()


@merge(buildBigWig,
       "bigwig_stats.load")
def loadBigWigStats(infiles, outfile):
    '''merge and load bigwig summary for all wiggle files.

    Summarise and merge bigwig files for all samples and load into a
    table called bigwig_stats

    Parameters
    ----------
    infiles : list
       Input filenames in :term:`bigwig` format
    outfile : string
        Output filename, the table name is derived from `outfile`.
    '''

    data = " ".join(
        ['<( bigWigInfo %s | perl -p -e "s/:/\\t/; s/ //g; s/,//g")' %
         x for x in infiles])
    headers = ",".join([P.snip(os.path.basename(x), ".bw")
                        for x in infiles])

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=track")

    statement = '''cgat combine_tables
    --header-names=%(headers)s
    --skip-titles
    --missing-value=0
    --ignore-empty
    %(data)s
    | perl -p -e "s/bin/track/"
    | cgat table2table --transpose
    | %(load_statement)s
    > %(outfile)s
    '''

    P.run()


@transform(MAPPINGTARGETS,
           regex(".bam"),
           ".bed.gz")
def buildBed(infile, outfile):
    ''' Generate :term:`bed` format file from :term:`bam` alignment file

    Parameters
    ----------
    infile : str
       Input filename in :term:`bam` format
    outfile : str
       Output filename in :term:`bed` format
    '''

    statement = '''
    cat %(infile)s
    | cgat bam2bed
          %(bed_options)s
          --log=%(outfile)s.log
          -
    | sort -k1,1 -k2,2n
    | bgzip
    > %(outfile)s;
    tabix -p bed %(outfile)s
    '''
    P.run()


@merge(buildBigWig, "igv_sample_information.tsv")
def buildIGVSampleInformation(infiles, outfile):
    '''build a file with IGV sample information

    Parameters
    ----------
    infiles : str
       Input filenames in :term:`bigwig` format
    outfile : str
       Output filename in :term:`tsv` format
    '''

    outf = IOTools.openFile(outfile, "w")
    first = True
    for fn in infiles:
        fn = os.path.basename(fn)
        parts = fn.split("-")
        if first:
            outf.write("sample\t%s\n" % "\t".join(
                ["%i" % x for x in range(len(parts))]))
            first = False
        outf.write("%s\t%s\n" % (fn, "\t".join(parts)))

    outf.close()


@follows(loadReadCounts,
         loadPicardStats,
         loadBAMStats,
         buildGeneProfiles,
         loadContextStats,
         loadIdxstats)
def general_qc():
    pass


@active_if(SPLICED_MAPPING)
@follows(loadExonValidation,
         loadGeneInformation,
         loadTranscriptLevelReadCounts,
         loadIntronLevelReadCounts,
         buildTranscriptProfiles,
         loadPicardRnaSeqMetrics,
         loadIdxstats)
def spliced_qc():
    pass


@follows(general_qc, spliced_qc)
def qc():
    pass


@follows(loadPicardDuplicationStats)
def duplication():
    pass


@follows(buildBed,
         loadBigWigStats,
         buildIGVSampleInformation)
def export():
    pass


@merge((loadBAMStats, loadPicardStats, loadContextStats), "view_mapping.load")
def createViewMapping(infile, outfile):
    '''create view in database for alignment stats.

    This view aggregates all information on a per-track basis.

    The table is built from the following tracks:

       context_stats
       bam_stats

    '''

    dbh = connect()

    tablename = P.toTable(outfile)
    view_type = "TABLE"

    tables = (("bam_stats", "track", ),
              ("context_stats", "track", ))

    # do not use: ("picard_stats_alignment_summary_metrics", "track"),)
    # as there are multiple rows per track for paired-ended data.

    P.createView(dbh, tables, tablename, outfile, view_type)


@follows(createViewMapping)
def views():
    pass


@follows(mapping, qc, views, duplication)
def full():
    pass


@follows(mapping)
def map_only():
    pass


@follows(mkdir("report"), mkdir(PARAMS.get("exportdir"), "export"))
def build_report():
    '''build report from scratch.'''
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''
    P.run_report(clean=False)


@follows(mkdir("%s/bamfiles" % PARAMS["web_dir"]),
         mkdir("%s/bigwigfiles" % PARAMS["web_dir"]),
         )
def publish():
    '''publish files.'''

    # directory, files
    export_files = {
        "bamfiles": glob.glob("*/*.bam") + glob.glob("*/*.bam.bai"),
        "bigwigfiles": glob.glob("*/*.bw"),
    }

    if PARAMS['ucsc_exclude']:
        for filetype, files in export_files.items():
            new_files = set(files)
            for f in files:
                for regex in P.asList(PARAMS['ucsc_exclude']):
                    if re.match(regex, f):
                        new_files.remove(f)
                        break

            export_files[filetype] = list(new_files)

    # publish web pages
    E.info("publishing report")
    P.publish_report(export_files=export_files)

    E.info("publishing UCSC data hub")
    P.publish_tracks(export_files)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
