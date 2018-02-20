"""
===================
Gtf_subset pipeline
===================

Overview
========

This pipeline generates a number of annotations that can be used with
downstream CGAT pipelines. The user will download a GTF from ENSEMBL
and then the GTF is parsed and filtered.In addition to downloading an
ensembl GTF the user will need to download an assembly report for their
specific genome and add it to the directory the pipeline is ran.

Common to all of the annotations generated in this pipeline is that they
are genomic - i.e. they are genomic intervals or relate to genomic intervals.
Thus, annotations are tied to a particular version of the genome. This is
parameterised within the pipeline.ini configuration file. The pipeline
follows two principle releases: the UCSC_ genome assembly and an ENSEMBL_
geneset version.

Note: This pipeline replaces pipeline_annotations which has now been moved
to the obsolete folder.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Principle targets
-----------------

lite
   This will run annotations that are used for our common upstream
   pipelines. i.e. pipeline_readqc.py, pipeline_mapping.py and
   pipeline_bamstats.py.

full
   This will run the entire pipeline


Configuration
-------------

The :file:`pipeline.ini` needs to be edited so that it points to the
appropriate locations of the auxiliary files.


On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Input
-----

This pipeline requires a Ensembl GTF, mirBase GFF3 file and an assembly report.

Ensembl GTF:
    This can be downloaded from http://www.ensembl.org/info/data/ftp/index.html
    Note: CGAT pipelines use the UCSC GTF convention (chr naming of contigs)
    and therefore the GTF is sanitized to the UCSC convention. As part of this
    process an ncbi assembly report needs to be specified (see below).

Assembly report:
    This is downloaded from the ncbi assembly page for your specific genome.
    Using hg19 as an example:
        Navigate to www......
        From the database tab select assembly and add your genome into the
        search bar i.e. hg19.
        Then click the link "Download the full sequence report"
        Add it to the folder where the pipeline will be ran, the file is
        for hg38 is called "GRCh38.p10_assembly_report.txt".

miRbase GFF3:
   This can be downloaded from miRbase http://www.mirbase.org/ftp.shtml.
   A path to the :term:`GFF3` file needs to be specified in the pipelin.ini
   configuration file. Make sure that the genome build version of the GFF3
   annotation file matches the ENSEMBL genome.

Running
-------
To run the pipeline, perform the following:

To run a basic set of annotations for the use of our common upstream pipelines
such as pipeline_readqc.py, pipeline_mapping.py and pipeline_bamstats.py you
can run the task `lite` as follows
::

   python /path/to/directory/pipeline_genesets.py make lite -v 5

To run the full set of annotations produced in this pipeline (which will be
used for some our our downstream pipelines such as pipeline_intervals.py)
you can run the `full` task:
::

   python /path/to/directory/pipeline_genesets.py make full -v 5


The pipeline can be run as any other CGAT pipeline, but as its purpose
is to provide a set of annotation that can be used by other pipelines
therefore there is an etiquette to be followed:

Using the pipeline results
--------------------------

The gtf_subset pipeline provides an interface for presenting its
results to other pipelines. The interface is defined in the file
:file:`pipeline.ini`. For example::

   [interface]
   # fasta file with cdna sequences
   cdna_fasta=ensembl.dir/cdna.fasta

The ini file of pipeline annotations can be loaded into the parameter
dictionary of your own pipeline::

    PARAMS.update(P.peekParameters(
         PARAMS["annotations_dir"],
         "pipeline_genesets.py",
         prefix="annotations_"),
         update_interface=True)

Parameters from the gtf_subset pipeline are now accessible via the
``annotations_`` prefix. As a result, the file
:file:`ensembl.dir/cdna.fasta` can be accessed as::

    PARAMS['annotations_cdna_fasta']


Working with non-ENSEMBL species
--------------------------------

:doc:`pipeline_gtf_subset` is very much wedded to annotations in ENSEMBL-
and UCSC_. Using a non-ENSEMBL species or non-UCSC species is possible by
building ENSEMBL- or UCSC-like input files. Even so, annotations that are
downloaded from the ENSEMBL or UCSC database will not be built. You will
thus need to ask if it is worth the effort.

As other pipelines will depend on the annotations in this pipeline it is
necessary to set up a :doc:`pipeline_gtf_subset` stub. To do so, simply
build the config files by running::

   python <SRC>pipeline_annotations.py config

and create the files that are being used in the downstream pipeline
explicitely (for example, for protein coding genes)::

   mkdir ensembl.dir
   cp <MYDATADIR>/my_gtf_geneset.gtf.gz ensembl.dir/geneset_coding.gtf.gz


Pipeline output
===============

The results of the computation are all stored in an sqlite relational
database file csvdb or as compressed files in genomic formats in the pipeline
directories. Output files are grouped by sections listed below.

The sections correspond to primary targets in the pipeline, i.e., to
build all annotations in the section ``assembly`` type::

   python <SRC>pipeline_annotations.py make assembly

Section: assembly
-----------------
Annotations derived from the genome assembly. Results are
in :file:`assembly.dir`.

contigs.tsv
   A :term:`tsv` formatted table with contig sizes

contigs.bed.gz
   bed file with contig sizes

section: ensembl
----------------

geneset_all.gtf.gz
   The full gene set after reconciling with assembly. Chromosomes names are
   renamed to be consistent with the assembly and some chromosomes
   are optionally removed. This file is the starting point for
   all annotations derived from the ENSEMBL geneset.

geneset_cds.gtf.gz
   A :term:`gtf` formatted file with only the CDS parts of transcripts.
   This set will naturally include only coding transcripts. UTR regions
   have been removed.

geneset_exons.gtf.gz
   A :term:`gtf` formatted file with only the exon parts of transcripts.
   This set includes both coding and non-coding transcripts. Coding
   transcripts span both the UTR and the CDS.

geneset_coding_exons.gtf.gz
   :term:`gtf` file with exon parts of protein coding transcripts.
   All other features are removed. These are all features annotated
   as "protein_coding" in the ENSEMBL gtf file.

geneset_noncoding_exons.gtf.gz
   :term:`gtf` file with exon parts of non-coding transcripts
   all other features are removed. These are all transcripts not
   annotated as "protein_coding" in the ENSEMBL gtf file.

geneset_lincrna_exons.gtf.gz
   :term:`gtf` file with exon parts of lincRNA transcripts. These
   are transcripts annotated as "lincRNA" in the ENSEMBL gtf file.

geneset_flat.gtf.gz
   A :term:`gtf` formatted file of flattened gene
   models. All overlapping transcripts have been merged. This set
   includes both coding and non-coding transcripts.

geneset_introns.gtf.gz
   A :term:`gtf` formatted file containing all intron features. All
   protein coding genes are retained and their exonic sequences are
   removed to retain introns from nested genes that may overlap.

section: mirbase
----------------

miRNA_non_primary_transcripts.gff3.gz
   A :term:`gff3` formatted file containing all of the non primary miRNA
   transcripts from mirbase

miRNA_primary_transcripts.gff3.gz
   A :term:`GFF3` formatted file containing all of the primery miRNA
   transcripts from miRbase.

section: ucsc
-------------

repeats.gff.gz
   :term:`gff` formatted file with structural/complex repeats

rna.gff.gz
   :term:`gff` formatted file with ribosomal rna annotations

section: geneset
----------------
Annotations derived from the ENSEMBL gene set. Annotations in
this section have been computed from the ENSEMBL gene set.
Results are in the directory :file:`geneset.dir`.

ref_flat.txt
   This creates a flat reference file from geneset_flat.gtf.gz
   for use in picard tools RNAseqmetrics.

section: bed
------------

This directory contains bed files that are generated from other annotations
in this pipeline.

genomic_context.bed.gz
   bed-formatted file with genomic context

Section: enrichment
-------------------

This section contains useful files for genomic enrichment analysis
a la gat_. The annotations are derived from other annotations in
this pipeline. Output files are in the directory :file:`enrichment.dir`.

annotation_gff.gz
   A :term:`gff` formatted file annotating the genome with respect
   to the geneset.  Annotations are non-overlapping and are based
   only on protein coding transcripts.

genestructure.gff.gz
   A :term:`gff` file annotation genomic regions by gene structure

territories.gff.gz
   gff file with gene territories, .i.e. regions around protein
   coding genes.  Intergenic space between genes is split at the
   midpoint between two genes.

tssterritories.gff.gz
   gff file with tss territories

greatdomains.gff.gz
   gff file of regulator domains defined a la GREAT

genomic_context_bed=genomic_context.bed.gz
   bed-formatted file with genomic context

genomic_function_bed=genomic_function.bed.gz
   bed-formatted file with functional annotations

genomic_function_tsv=genomic_function.tsv.gz
   tsv-formatted file mapping terms to descriptions


Database design
---------------

Tables in the database usually represent genomic features such as
transcripts, genes or chromosomes. These are identified by the
following columns:

+--------------------+-----------------------------------------+
|*Column*            |*Content*                                |
+--------------------+-----------------------------------------+
|transcript_id       |ENSEMBL transcript identifier            |
+--------------------+-----------------------------------------+
|gene_id             |ENSEMBL gene id                          |
+--------------------+-----------------------------------------+
|contig              |Chromosome name                          |
+--------------------+-----------------------------------------+

Example
=======

**Supply example data**


====
Code
====
"""
import sys
import re
import os
import sqlite3
import glob
import pandas as pd
from ruffus import follows, transform, merge, mkdir, files, jobs_limit,\
    suffix, regex, add_inputs, originate
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineGtfsubset as PipelineGtfsubset
import CGATPipelines.PipelineUCSC as PipelineUCSC
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineGO as PipelineGO

###################################################
# Pipeline configuration
###################################################
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# Add automatically created files to the interface.  This is required
# when the pipeline is peek'ed.  The statement below will
# add the following to the dictionary:
#
# "geneset.dir/lincrna_gene_tss.bed.gz" maps to
# "interface_geneset_lincrna_gene_tss_bed"
PARAMS.update(dict([
    ("interface_geneset_%s" %
     re.sub("[.]", "_", os.path.basename(P.snip(x, ".gz"))), x)
    for x in glob.glob('geneset.dir/*.bed.gz')]))


def connect():
    '''connect to database.'''

    dbh = sqlite3.connect(PARAMS["database_name"])
    return dbh


def connectToUCSC():
    return PipelineGtfsubset.connectToUCSC(
        host=PARAMS["ucsc_host"],
        user=PARAMS["ucsc_user"],
        database=PARAMS["ucsc_database"])

############################################################
# Assembly
############################################################


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_contigs'])
def buildContigSizes(infile, outfile):
    '''
    Get contig sizes from indexed genome :term:`fasta` files and
    outputs to a text file.
    Parameters
    ----------
    infile : str
      infile is constructed from the `PARAMS` variable to retrieve
      the `genome` :term:`fasta` file
    Returns
    -------
    outfile : str
      outfile is a text format file that contains two columns, matched
      contig name and contig size (in nucleotides).  The output file
      name is defined in `PARAMS: interface_contigs`.
    '''

    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    contigs = []

    for contig, size in fasta.getContigSizes(with_synonyms=False).items():
        contigs.append([contig, size])
    df_contig = pd.DataFrame(contigs, columns=['contigs', 'size'])
    df_contig.sort_values('contigs', inplace=True)
    df_contig.to_csv(outfile, sep="\t", header=False, index=False)


@follows(buildContigSizes)
@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_contigs_bed'])
def buildContigBed(infile, outfile):
    '''
    Gets the contig sizes and co-ordinates from an indexed genome :term:`fasta`
    file and outputs them to :term:`BED` format
    Parameters
    ----------
    infile : str
      infile is constructed from `PARAMS` variable to retrieve
      the `genome` :term:`fasta` file
    Returns
    -------
    outfile : str
      :term:`BED` format file containing contig name, value (0) and contig size
      in nucleotides.  The output file name is defined in
      `PARAMS: interface_contigs_bed`
    '''
    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs = IOTools.openFile(outfile, "w")

    for contig, size in fasta.getContigSizes(with_synonyms=False).items():
        outs.write("%s\t%i\t%i\n" % (contig, 0, size))

    outs.close()


@follows(buildContigBed)
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       (PARAMS['interface_contigs_ungapped_bed'],
        PARAMS['interface_gaps_bed'],
        ))
def buildUngappedContigBed(infile, outfiles):
    '''
    Constructs :term:`BED` format files containing both gapped and ungapped
    contig sizes from an index genome :term:`fasta` file.

    Parameters
    ----------
    infile: str
      infile is constructed from `PARAMS` variable to retrieve
      the `genome` :term:`fasta` file

    assembly_gaps_min_size: int
      `PARAMS` - the minimum size (in nucleotides) for an assembly gap

    Returns
    -------
    outfiles: list
      two separate :term:`BED` format output files containing the contig sizes
      for contigs with and without gaps.  The names are defined
      in the `PARAMS` `interface_contigs_ungapped_bed` and
      `interface_gaps_bed` parameters.
    '''

    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs_nogap = IOTools.openFile(outfiles[0], "w")
    outs_gap = IOTools.openFile(outfiles[1], "w")
    min_gap_size = PARAMS["assembly_gaps_min_size"]

    for contig, size in fasta.getContigSizes(with_synonyms=False).items():

        seq = fasta.getSequence(contig)

        def gapped_regions(seq):
            is_gap = seq[0] == "N"
            last = 0
            for x, c in enumerate(seq):
                if c == "N":
                    if not is_gap:
                        last = x
                        is_gap = True
                else:
                    if is_gap:
                        yield(last, x)
                        last = x
                        is_gap = False
            if is_gap:
                yield last, size

        last_end = 0
        for start, end in gapped_regions(seq):
            if end - start < min_gap_size:
                continue

            if last_end != 0:
                outs_nogap.write("%s\t%i\t%i\n" % (contig, last_end, start))
            outs_gap.write("%s\t%i\t%i\n" % (contig, start, end))
            last_end = end

        if last_end < size:
            outs_nogap.write("%s\t%i\t%i\n" % (contig, last_end, size))

    outs_nogap.close()
    outs_gap.close()


@follows(buildUngappedContigBed)
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_cpg_bed'])
def buildCpGBed(infile, outfile):
    '''
    Output a :term:`BED` file that contains the location of all CpGs
    in the input genome using `CGAT` script `fasta2bed`.

    Parameters
    ----------
    infile: str
      infile is constructed from `PARAMS` variable to retrieve
      the `genome` :term:`fasta` file

    Returns
    -------
    outfile: str
      A :term:`BED` format file containing location of CpGs across the
      genome.  The BED file is then indexed using tabix
    '''

    job_memory = PARAMS["job_highmemory"]

    statement = '''
    cgat fasta2bed
        --method=cpg
        --log=%(outfile)s.log
    < %(infile)s
    | bgzip
    > %(outfile)s
    '''

    P.run()

    statement = '''
    tabix -p bed %(outfile)s
    '''
    P.run()

###################################################################
# ENSEMBL gene set
###################################################################


@follows(mkdir('ensembl.dir'))
@transform(PARAMS["ensembl_filename_gtf"],
           regex("(\S+)"),
           r"%s" % PARAMS['interface_geneset_all_gtf'])
def buildUCSCGeneSet(infile, outfile):
    '''output sanitized ENSEMBL geneset.

    This method outputs an ENSEMBL gene set after some sanitizing steps:

    1. Chromosome names are changed to the UCSC convention.
    2. Chromosomes that match the regular expression specified in
       the configuration file are removed.

    Arguments
    ---------
    infiles : string
       ENSEMBL geneset in :term:`gtf` format.
       NCBI Assembly report in 'txt' format.
    outfile : string
       geneset in :term:`gtf` format.

    '''

    job_memory = PARAMS["job_memory"]

    statement = ['''zcat %(infile)s
    | grep 'transcript_id'
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=ucsc
    --assembly-report="%(ncbi_assembly_report)s"
    --log=%(outfile)s.log
    ''']

    if PARAMS["ncbi_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --contig-pattern="%(ncbi_remove_contigs)s" ''')

    statement.append(
        '''
        | cgat gtf2gtf
        --method=set-gene_biotype-to-source
        --log=%(outfile)s.log
        | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_cds_gtf'])
def buildCdsTranscript(infile, outfile):
    '''
    Output the CDS features from an ENSEMBL gene set

    takes all of the features from a :term:`gtf` file
    that are feature types of ``CDS``.

    Note - we have not filtered on gene_biotype because some of the CDS
    are classified as polymorphic_pseudogene.

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''

    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["CDS"]

    m.filterGTF(outfile, filteroption, filteritem, operators=None)


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_exons_gtf'])
def buildExonTranscript(infile, outfile):
    '''
    Output of the exon features from an ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``exon``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["exon"]

    m.filterGTF(outfile, filteroption, filteritem, operators=None)


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_coding_exons_gtf'])
def buildCodingExonTranscript(infile, outfile):
    '''
    Output of the coding exon features from abn ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``exon``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroption = [PARAMS['ensembl_cgat_feature'],
                    PARAMS['ensembl_cgat_gene_biotype']]
    filteritem = ["exon", "protein_coding"]

    m.filterGTF(outfile, filteroption, filteritem, operators="and")


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_lincrna_exons_gtf'])
def buildLincRNAExonTranscript(infile, outfile):
    '''
    Output of the lincRNA features from an ENSEMBL gene set

    Takes all of the features from a :term:`gtf` file
    that are features of ``lincRNA``

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroptions = [PARAMS['ensembl_cgat_feature'],
                     PARAMS['ensembl_cgat_gene_biotype']]

    filteritem = ["exon", "lincRNA"]

    m.filterGTF(outfile, filteroptions, filteritem, operators="and")


@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_noncoding_exons_gtf'])
def buildNonCodingExonTranscript(infile, outfile):
    '''
    Output of the non-coding exon features from an ENSEMBL gene set

    Remove all of the features from a :term:`gtf` file
    that are features of ``exon`` and are protein-coding

    Arguments
    ---------
    infile : from ruffus
       ENSEMBL geneset, filename named in pipeline.ini
    outfile : from ruffus
       Output filename named in pipeline.ini
    filteroption : string
       Filter option set in the piepline.ini as feature column in GTF
       nomenclature
    '''
    m = PipelineGtfsubset.SubsetGTF(infile)

    filteroptions = [PARAMS['ensembl_cgat_feature'],
                     PARAMS['ensembl_cgat_gene_biotype']]
    filteritem = ["exon", "protein_coding"]

    m.filterGTF(outfile, filteroptions, filteritem, operators="and not")


@transform((buildUCSCGeneSet,
            buildCdsTranscript,
            buildExonTranscript,
            buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           suffix(".gtf.gz"), "_gtf.load")
def loadTranscripts(infile, outfile):
    '''load transcripts from a GTF file into the database.

    The table will be indexed on ``gene_id`` and ``transcript_id``

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''

    job_memory = PARAMS["job_highmemory"]

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id "
        "--allow-empty-file ")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2tsv -f
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


@P.add_doc(PipelineGtfsubset.buildFlatGeneSet)
@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS['interface_geneset_flat_gtf'])
def buildFlatGeneSet(infile, outfile):
    PipelineGtfsubset.buildFlatGeneSet(infile, outfile,
                                       job_memory=PARAMS["job_highmemory"])

####################################################################
# Geneset derived annotations
####################################################################


@follows(mkdir("geneset.dir"))
@transform(buildUCSCGeneSet,
           suffix("ensembl.dir/geneset_all.gtf.gz"),
           PARAMS["interface_ref_flat"])
def buildRefFlat(infile, outfile):
    '''build flat geneset for Picard RnaSeqMetrics.
    '''

    tmpflat = P.getTempFilename(".")

    job_memory = PARAMS["job_memory"]

    statement = '''
    gtfToGenePred -genePredExt -geneNameAsName2 %(infile)s %(tmpflat)s;
    paste <(cut -f 12 %(tmpflat)s) <(cut -f 1-10 %(tmpflat)s)
    > %(outfile)s
    '''
    P.run()
    os.unlink(tmpflat)


@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transcript_region.bed.gz')
def buildTranscriptRegions(infile, outfile):
    """
    export a table of seleno cysteine transcripts.

    Selenocysteine containing transcripts are identified by checking
    if their protein sequence contains ``U``.

    The table contains a single column ``transcript_id`` with ENSEMBL
    transcript identifiers as values.
    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """
    gunzip < %(infile)s
    | cgat gtf2gtf --method=join-exons
    --log=%(outfile)s.log
    | cgat gff2bed --is-gtf
    --set-name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_region.bed.gz')
def buildGeneRegions(infile, outfile):
    """build a :term:`bed` file of regions spanning whole gene models.

    This method outputs a single interval spanning the genomic region
    that covers all transcripts within a particular gene.

    The name column of the :term:`bed` file is set to the `gene_id`.

    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=merge-transcripts
    --log=%(outfile)s.log
    | cgat gff2bed --is-gtf --set-name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@follows(mkdir("geneset.dir"))
@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transcript_tss.bed.gz')
def buildTranscriptTSS(infile, outfile):
    """build a :term:`bed` file with transcription start sites.

    This method outputs all transcription start sites within a
    geneset. The trancription start site is derived from the most
    upstream coordinate of each transcript.

    The name column of the :term:`bed` file is set to the
    `transcript_id`.

    Arguments
    ---------
    infile : list
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """
    gunzip < %(infile)s
    | cgat gtf2gtf --method=join-exons
    --log=%(outfile)s.log
    | cgat gtf2gff --method=promotors
    --promotor-size=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | cgat gff2bed --is-gtf --set-name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tssinterval.bed.gz')
def buildGeneTSSInterval(infile, outfile):
    """build a :term:`bed` file with intervals that cover all transcription
    start sites within a gene.

    This method outputs for each gene the smallest genomic region that covers
    all the transcription start sites within that gene.

    The name column of the :term:`bed` file is set to the
    `gene_id`.

    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=join-exons
    --log=%(outfile)s.log
    | cgat gtf2gff
    --method=promotors
    --promotor-size=1
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | sed s/transcript/exon/g
    | sed s/exon_id/transcript_id/g
    | cgat gtf2gtf
    --method=merge-transcripts
    --log=%(outfile)s.log
    | cgat gff2bed
    --is-gtf
    --set-name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transcript_tts.bed.gz')
def buildTranscriptTTS(infile, outfile):
    """build a :term:`bed` file with transcription termination sites.

    This method outputs all transcription start sites within a
    geneset. The trancription start site is derived from the most
    downstream coordinate of each transcript.

    The name column of the :term:`bed` file is set to the
    `transcript_id`.

    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """
    gunzip < %(infile)s
    | cgat gtf2gtf --method=join-exons
    --log=%(outfile)s.log
    | cgat gtf2gff --method=tts
    --promotor-size=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | cgat gff2bed --is-gtf --set-name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@follows(mkdir("geneset.dir"))
@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tss.bed.gz')
def buildGeneTSS(infile, outfile):
    """build a :term:`bed` file with transcription start sites per gene.

    This method outputs a single transcription start sites for each
    gene within a geneset. The trancription start site is derived from
    the most upstream coordinate of each gene.

    The name column of the :term:`bed` file is set to the
    `gene_id`.

    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """gunzip < %(infile)s
    | cgat gtf2gtf
    --method=merge-transcripts
    --log=%(outfile)s.log
    | cgat gtf2gff --method=promotors --promotor-size=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | cgat gff2bed --is-gtf --set-name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s"""
    P.run()


@transform((buildCodingExonTranscript,
            buildNonCodingExonTranscript,
            buildLincRNAExonTranscript),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tts.bed.gz')
def buildGeneTTS(infile, outfile):
    """build a :term:`bed` file with transcription termination sites per gene.

    This method outputs a single transcription start sites for each
    gene within a geneset. The trancription start site is derived from
    the most downstream coordinate of each gene.

    The name column of the :term:`bed` file is set to the
    `gene_id`.

    Arguments
    ---------
    infile : string
       Input filename with geneset in :term:`gtf` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.

    """

    job_memory = PARAMS["job_memory"]

    statement = """gunzip < %(infile)s
    | cgat gtf2gtf
    --method=merge-transcripts
    --log=%(outfile)s.log
    | cgat gtf2gff --method=tts --promotor-size=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | cgat gff2bed --is-gtf --set-name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s"""
    P.run()


@transform(buildGeneRegions,
           regex('(.*)_.*.bed.gz'),
           add_inputs(buildContigSizes),
           r'\1_intergenic.bed.gz')
def buildIntergenicRegions(infiles, outfile):
    """build a :term:`bed` file with regions not overlapping any genes.

    Arguments
    ---------
    infiles : list
       - Input filename with geneset in :term:`gtf` format.
       - Input filename with chromosome sizes in :term:`tsv` format.
    outfile : string
       Output filename with genomic regions in :term:`bed` format.
    """

    infile, contigs = infiles

    job_memory = PARAMS["job_memory"]

    statement = '''zcat %(infile)s
    | sort -k1,1 -k2,2n
    | complementBed -i stdin -g %(contigs)s
    | gzip
    > %(outfile)s'''
    P.run()


@P.add_doc(PipelineGtfsubset.loadGeneInformation)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(mkdir('ensembl.dir'))
@transform(PARAMS["ensembl_filename_gtf"],
           suffix(PARAMS["ensembl_filename_gtf"]),
           "ensembl.dir/gene_info.load")
def loadGeneInformation(infile, outfile):
    '''load the transcript set.'''
    PipelineGtfsubset.loadGeneInformation(infile, outfile,
                                          job_memory=PARAMS["job_highmemory"])


@follows(loadGeneInformation)
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

    select = dbh.execute("""SELECT DISTINCT gene_id
    FROM gene_info
    WHERE gene_biotype = 'protein_coding'""" % locals())

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("gene_id\n")
        outf.write("\n".join((x[0] for x in select)) + "\n")


@transform(buildUCSCGeneSet,
           regex(".*"),
           PARAMS['interface_utr_all_gtf'])
def buildUtrGeneSet(infile, outfile):

    job_memory = PARAMS["job_memory"]

    statement = "zcat %(infile)s | grep 'utr' | gzip > %(outfile)s"

    P.run()


@transform(buildFlatGeneSet,
           regex(".*"),
           add_inputs(identifyProteinCodingGenes,
                      buildExonTranscript),
           PARAMS['interface_geneset_intron_gtf'])
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

    infile, genes_tsv, filename_exons = infiles

    job_memory = PARAMS["job_memory"]

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
    | awk -v OFS="\\t" -v FS="\\t" '{$3="intron"; print}'
    | gzip
    > %(outfile)s
    '''
    P.run()

# Next need to add identifyProteinCodingGenes, buildIntronGeneModels
# aim is to generate the intron gtf here for use in bamstats

################################################################
# UCSC derived annotations
################################################################


@follows(mkdir('ucsc.dir'))
@originate(PARAMS["interface_rna_gff"])
def importRNAAnnotationFromUCSC(outfile):
    """This task downloads UCSC repetetive RNA types.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_rnatypes"]),
        outfile=outfile,
        remove_contigs_regex=PARAMS["ncbi_remove_contigs"],
        job_memory=PARAMS["job_memory"])


@follows(mkdir('ucsc.dir'))
@originate(PARAMS["interface_repeats_gff"])
def importRepeatsFromUCSC(outfile):
    """This task downloads UCSC repeats types as identified
    in the configuration file.
    """
    PipelineGtfsubset.getRepeatDataFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_repeattypes"]),
        outfile=outfile,
        job_memory=PARAMS["job_memory"])


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((importRepeatsFromUCSC,
            importRNAAnnotationFromUCSC),
           suffix(".gff.gz"), "_gff.load")
def loadRepeats(infile, outfile):
    """load genomic locations of repeats into database.

    This method loads the genomic coordinates (contig, start, end)
    and the repeat name into the database.

    Arguments
    ---------
    infile : string
        Input filename in :term:`gff` with repeat annotations.
    outfile : string
        Output filename with logging information. The table name is
        derived from outfile.

    """

    job_memory = PARAMS["job_memory"]

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=class "
        "--header-names=contig,start,stop,class")

    statement = """zcat %(infile)s
    | cgat gff2bed --set-name=class
    | grep -v "#"
    | cut -f1,2,3,4
    | %(load_statement)s
    > %(outfile)s"""
    P.run()


#######################################################
# miRBase annotations
########################################################

@follows(mkdir("mirbase.dir"))
@transform(PARAMS['mirbase_filename_mir_gff'],
           suffix(PARAMS['mirbase_filename_mir_gff']),
           PARAMS['interface_geneset_primary_mir_gff'])
def buildmiRPrimaryTranscript(infile, outfile):

    '''
    This function will subset a miRbase annotation gff3 file.The GFF3
    file can be downloaded from miRbase. Make sure the annotation matches
    the genome build that you are using.

    This function will subset the GFF3 file by selecting annotations that are
    labled "miRNA_primary_transcript"
    '''

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA_primary_transcript"]

    m.filterGFF3(outfile, filteroption, filteritem)


@follows(buildmiRPrimaryTranscript)
@transform(PARAMS['mirbase_filename_mir_gff'],
           suffix(PARAMS['mirbase_filename_mir_gff']),
           PARAMS['interface_geneset_mir_gff'])
def buildmiRNonPrimaryTranscript(infile, outfile):

    '''
    This function will subset a miRbase annotation gff3 file.The GFF3
    file can be downloaded from miRbase. Make sure the annotation matches
    the genome build that you are using.

    This function will subset the GFF3 file by selecting annotations that are
    labled "miRNA". This will subset all of the non primary transcripts.
    '''

    m = PipelineGtfsubset.SubsetGFF3(infile)

    filteroption = PARAMS['ensembl_cgat_feature']
    filteritem = ["miRNA"]

    m.filterGFF3(outfile, filteroption, filteritem)


@transform((buildmiRPrimaryTranscript,
            buildmiRNonPrimaryTranscript),
           suffix(".gff3.gz"), "_gff3.load")
def loadmiRNATranscripts(infile, outfile):
    '''load transcripts from a GFF3 file into the database.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gff3` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''

    job_memory = PARAMS["job_memory"]

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--allow-empty-file "
        "--header-names=feature,Name")

    statement = '''
     export LANG=en_GB.UTF-8 && zcat %(infile)s
    | cgat gtf2tsv --is-gff3 --attributes-as-columns 2> /dev/null
    | grep -v "#"
    | cut -f3,12
    |%(load_statement)s
    > %(outfile)s'''
    P.run()

###############################################################
# Ontologies
###############################################################


@P.add_doc(PipelineGO.createGOFromENSEMBL)
@follows(mkdir('ontologies.dir'))
@files([(None, PARAMS["interface_go_ensembl"]), ])
def createGO(infile, outfile):
    '''
    Downloads GO annotations from ensembl
    Uses the go_host, go_database and go_port parameters from the ini file
    and runs the runGO.py "filename-dump" option.
    This calls DumpGOFromDatabase from GO.py
    '''
    PipelineGO.createGOFromENSEMBL(infile, outfile,
                                   job_memory=PARAMS["job_highmemory"])


@P.add_doc(PipelineGO.createGOSlimFromENSEMBL)
@transform(createGO,
           regex("(.*)"),
           PARAMS["interface_goslim_ensembl"])
def createGOSlim(infile, outfile):
    '''
    Downloads GO slim annotations from ensembl
    '''
    E.warn(PARAMS['go_url_goslim'])
    PipelineGO.createGOSlimFromENSEMBL(infile, outfile,
                                       job_memory=PARAMS["job_highmemory"])


@P.add_doc(P.load)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((createGO, createGOSlim),
           suffix(".tsv.gz"),
           r"\1_assignments.load")
def loadGOAssignments(infile, outfile):
    '''
    Load GO assignments into database.'''
    P.load(infile, outfile,
           options="--add-index=gene_id --add-index=go_id")


################################################################
# Enrichment analysis
#################################################################

@P.add_doc(PipelineGeneset.annotateGenome)
@follows(mkdir('enrichment.dir'))
@files(buildUCSCGeneSet, PARAMS['interface_annotation_gff'])
def annotateGenome(infile, outfile):
    """This task only considers protein coding genes as
    processed_transcripts tend to cover larger genomic regions and
    often overlap between adjacent protein coding genes.

    """
    PipelineGeneset.annotateGenome(infile,
                                   outfile,
                                   only_proteincoding=True,
                                   job_memory=PARAMS["job_memory"])


@P.add_doc(PipelineGeneset.annotateGeneStructure)
@follows(mkdir('enrichment.dir'))
@files(buildUCSCGeneSet, PARAMS['interface_genestructure_gff'])
def annotateGeneStructure(infile, outfile):
    """This task only considers protein coding genes as
    processed_transcripts tend to cover larger genomic regions and
    often overlap between adjacent protein coding genes.

    """
    PipelineGeneset.annotateGeneStructure(infile,
                                          outfile,
                                          only_proteincoding=True,
                                          job_memory=PARAMS["job_memory"])


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_territories_gff"])
def buildGeneTerritories(infile, outfile):
    """build gene territories from protein coding genes.

    The territory of a gene is defined as the region of the
    gene extended by a certain radius on either end. If the
    gene territories of two genes overlap, they are resolved
    at the mid-point between the two adjacent genes.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gff` format.
    enrichment_territories_radius : int
       see :term:`PARAMS`
    """

    job_memory = PARAMS["job_highmemory"]

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter
    --filter-method=proteincoding
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort --sort-order=gene
    | cgat gtf2gtf
    --method=merge-transcripts
    | cgat gtf2gtf
    --method=sort --sort-order=position
    | cgat gtf2gff
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --territory-extension=%(enrichment_territories_radius)s
    --method=territories
    | cgat gtf2gtf
    --method=filter
    --filter-method=longest-gene
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s '''

    P.run()


@P.add_doc(PipelineGeneset.buildGenomicFunctionalAnnotation)
@follows(mkdir('enrichment.dir'))
@merge((buildGeneTerritories, loadGOAssignments),
       (PARAMS["interface_genomic_function_bed"],
        PARAMS["interface_genomic_function_tsv"],
        ))
def buildGenomicFunctionalAnnotation(infiles, outfiles):

    territories_gtf_file = infiles[0]

    PipelineGeneset.buildGenomicFunctionalAnnotation(
        territories_gtf_file,
        dbh=connect(),
        outfiles=outfiles,
        job_memory=PARAMS["job_memory"])


@P.add_doc(PipelineGtfsubset.buildGenomicContext)
@follows(mkdir('enrichment.dir'))
@merge((importRepeatsFromUCSC,
        importRNAAnnotationFromUCSC,
        buildUCSCGeneSet,
        buildUtrGeneSet,
        buildIntronGeneModels),
       PARAMS["interface_genomic_context_bed"])
def buildGenomicContext(infiles, outfile):
    PipelineGtfsubset.buildGenomicContext(infiles, outfile,
                                          job_memory=PARAMS["job_highmemory"])


##############################################################
# Define final task
###############################################################

@follows(buildUCSCGeneSet,
         buildContigSizes,
         buildExonTranscript,
         buildCodingExonTranscript,
         buildFlatGeneSet,
         buildIntronGeneModels,
         buildGenomicContext,
         buildRefFlat,
         importRNAAnnotationFromUCSC)
def lite():
    '''
    build only tasks that are used by the common upstream pipeline,
    pipeline_readqc, mapping and bamstats. To run the downstream pipelines
    run the full task so all annotations are made.
    '''
    pass


@follows(buildGenomicContext,
         buildGenomicFunctionalAnnotation)
def enrichment():
    """convenience target : annotations for enrichment analysis"""


@follows(buildTranscriptRegions,
         buildTranscriptTSS,
         buildTranscriptTTS,
         buildGeneRegions,
         buildGeneTSS,
         buildGeneTTS,
         buildGeneTSSInterval,
         buildIntergenicRegions)
def geneset():
    """convenience target : geneset derived annotations"""


@follows(loadTranscripts,
         loadRepeats,
         loadmiRNATranscripts,
         loadGeneInformation,
         buildFlatGeneSet,
         buildRefFlat,
         buildUtrGeneSet,
         buildContigBed,
         buildGeneTerritories,
         geneset,
         enrichment)
def full():
    '''build all targets - A dummy task to run the pipeline to
    completion.'''
    pass


@follows(mkdir("Report.dir"))
def build_report():
    '''report dummy task'''
    pass


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
