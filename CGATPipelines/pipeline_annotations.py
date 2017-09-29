"""===================
Annotation pipeline
===================

The annotation pipeline imports various third party annotations
or creates them for use in other pipelines.

The purpose of this pipeline is to automate and standardize the
way we retrieve and build genomic annotations but also to allow
sharing of annotations between projects and people. An important
part is the reconciliation of different data sources in terms
of chromosome names.

Common to all annotations in this pipeline is that they are genomic -
i.e. they are genomic intervals or relate to genomic intervals. Thus,
annotations are tied to a particular version of a genome. This pipeline
follows two principal releases: the UCSC_ genome assembly version and an
ENSEMBL_ geneset version.

The pipeline contains multiple sections that can be built on demand
or when relevant. Certain annotations (ENCODE, GWAS data) exist only
for specific species. The sections are:

assembly
   Genome assembly related information such as the location of
   gaps, chromosome lengths, etc.

ucsc
   Typical annotations downloaded from UCSC such as repeats.

ensembl
   The Ensembl gene set, reconciled with the assembly,
   and various subsets (coding genes, noncoding genes, ...).

geneset
   Annotations derived from the ENSEMBL gene set.

enrichment
   Annotations of genomic regions useful for enrichment
   analysis. These are derived from multiple input sources.

gwas
   GWAS data from the GWAS Catalog and DistlD

ontologies
   Ontology annotations (GO, KEGG) of genes.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The :file:`pipeline.ini` needs to be edited so that it points to the
appropriate locations of the auxiliary files. See especially:

1 section ``[ensembl]`` with the location of the ENSEMBL dump
    files (``filename_gtf``, filename_pep``, ``filename_cdna``)

2 section ``[general]`` with the location of the indexed genomic
    fasta files to use and the name of the genome, as well as the
    genome assembly report obtained from NCBI for mapping between
    UCSC and ENSEMBL contigs. This can be obtained from:
    https://www.ncbi.nlm.nih.gov/assembly
    see :doc:`../modules/IndexedFasta`.

3 section ``[ucsc]`` with the name of the database to use (default=``hg19``).

Input
-----

This script requires no input within the :term:`working directory`, but
will look up some files in directories specified in the configuration
file :file:`pipeline.ini` and download annotations using mysql.

Running
-------

The pipeline can be run as any other CGAT pipeline, but as its purpose
is to provide a set of shared annotation between multiple projects
there is an etiquette to be followed:

Using the pipeline results
--------------------------

The annotations pipeline provides an interface for presenting its
results to other pipelines. The interface is defined in the file
:file:`pipeline.ini`. For example::

   [interface]
   # fasta file with cdna sequences
   cdna_fasta=ensembl.dir/cdna.fasta

The ini file of pipeline annotations can be loaded into the parameter
dictionary of your own pipeline::

    PARAMS.update(P.peekParameters(
         PARAMS["annotations_dir"],
         "pipeline_annotations.py",
         prefix="annotations_"),
         update_interface=True)

Parameters from the annotation pipeline are now accessible via the
``annotations_`` prefix. As a result, the file
:file:`ensembl.dir/cdna.fasta` can be accessed as::

    PARAMS['annotations_cdna_fasta']

Extending the pipeline
-----------------------

Please feel free to add more annotations to the pipeline, but
considering its shared usage, please consult with others. In
particular, consider the following questions:

1. Is the annotation that I want to add genomic? For example,
   protein-protein interaction data should be organized separately.

2. Is the annotation of general interest? Do not add if an annotation
   is specific to a particular species or of very specialized
   interest. Note that there are some exceptions for annotations from
   certain species (human).

3. Is the annotation subjective? The pipeline consciously
   avoids providing annotations for regions such as promotors as their
   definition varies from person to person. Instead, the pipeline
   presents files with unambiguous coordinates such as transcription
   start sites. In the case of promotors, these could be derived from
   transcription start sites and the ``bedtools extend`` command.

4. What is the right format for the annotation? :term:`bed` formatted
   file are ideal for intervals with a single annotation. If multiple
   annotations are assigned with a feature, use :term:`gff`. For genes,
   use :term:`gtf`. Do not provide the same information with different
   formats - formats can be easily interconverted using CGAT tools.

Known problems
--------------

The pipeline takes its basic information about the genome and genes
from files downloaded from the genome browsers:

* UCSC: the genomic sequence in :term:`fasta` format.
* ENSEMBL: the gene set in :term:`GTF` format.

Additional data is downloaded from the genome browser databases either
via mysql or through biomart. It is thus important that the releases of
this additional data is consistent with the input files above.

.. note::

    The mechanism for getting biomart to download data for
    a particular ENSEMBL release involves changing the biomart server
    to an archive server.

Also, other data sources will have release cycles that are not tied
to a particular UCSC or ENSEMBL release. It is important to coordinate
and check when updating these other data sources.

Working with non-ENSEMBL species
--------------------------------

:doc:`pipeline_annotations` is very much wedded to annotations in ENSEMBL-
and UCSC_. Using a non-ENSEMBL species or non-UCSC species is possible by
building ENSEMBL- or UCSC-like input files. Even so, annotations that are
downloaded from the ENSEMBL or UCSC database will not be built. You will
thus need to ask if it is worth the effort.

As many other pipelines depend on the annotations in this pipeline it is
necessary to set up a :doc:`pipeline_annotations` stub. To do so, simply
build the config files by running::

   python <SRC>pipeline_annotations.py config

and create the files that are being used in the downstream pipeline
explicitely (for example, for protein coding genes)::

   mkdir ensembl.dir
   cp <MYDATADIR>/my_gtf_geneset.gtf.gz ensembl.dir/geneset_coding.gtf.gz

Roadmap
-------

There are many annotations that could possibly be brought into this pipeline:

* ENCODE data
     Can be used directly from a download directory?

* Genome segmentation based on ENCODE
     Definitions of enhancers, etc. Note that these will depend not on the
     genome, but on the cell type as well and thus might be project specific?

* Gene networks
     Functional assocation between genes. Outside of the
     scope of this pipeline?

* Mapability
     Mapability tracks are not available from all genomes. The pipeline
     could include runnig GEM on the assembly. For now it has been taken
     out as it is a rather long job.

Pipeline output
===============

The results of the computation are all stored in an sqlite relational
database file or as compressed files in genomic formats in the pipeline
directory. Output files are grouped by sections listed below.

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

contigs_ungapped.bed.gz
   :term:`bed` file with contigs excluding any gapped regions

gaps.bed.gz
   :term:`bed` file with gapped regions in contigs

genome.tsv.gz
   chromosome nucleotide composition and other stats

cpg.bed.gz
   filename with locations of CpG in bed format

gc_segmentation.bed.gz
   bed file with genome segmented into regions of similar G+C content
   using naive window based classification

gcprofile_bins.bed.gz
   bed file with genome segmented according to similar G+C content
   using the GCProfile method

Tables:
   genome
      Nucleotide composition of chromosomes

Section: ucsc
-------------

Various UCSC derived annotations. Results are in the
:file:`ucsc.dir`.

cpgislands.bed.gz
    :term:`bed` file with gapped regions in contigs

repeats.gff.gz
    :term:`gff` formatted file with structural/complex repeats

allrepeats.gff.gz
    :term:`gff` formatted file with all repeats including
    simple repeats

rna.gff.gz
    :term:`gff` formatted file with ribosomal rna annotations

repeats.gff.gz
    A :term:`gff` formatted file of repetitive sequences (obtained
    from UCSC repeatmasker tracks).

rna.gff.gz
    A :term:`gff` formatted file of repetitive RNA sequences in the genome
    (obtained from UCSC repeatmasker tracks).

mapability_xx.bed.gz
    Mapability files from UCSC CRG Alignability tracks. XX is the read
    length.

mapability_xx.bed.filtered.gz
    Similar to mapability_xx.bed.gz, but short regions of low mapability
    have been merged.

Tables
   repeat
       complex repeat locations
   repeat_counts
       Number of occurances for each repeat type.


Section: ensembl
----------------

Annotations within the ENSEMBL gene set after reconciliation
with the UCSC genome assembly. The results are in :file:`ensembl.dir`.
Annotations here are the original ENSEMBL annotations bar some
filtering.

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

peptides.fasta
   A :term:`fasta` formatted file of peptide sequences of coding
   transcripts.

cds.fasta
   A :term:`fasta` formatted file of coding sequence of coding transcripts.

cdna.fasta
   A :term:`fasta` formatted file of transcripts including both
   coding and non-coding parts.

Tables:
   transcript_info
       Information about transcripts (gene, biotype, status, ...)
       downloaded from biomart.

   transcript_synonyms
       Alternative names for transcripts

   gene_info
       Information about ENSEMBL genes in ENSEMBL gtf file.

   ensembl_to_entrez
       Table mapping ENSEMBL gene identifiers to
       ENTREZ identifiers

   cds_stats
       Table with nucleotide composition of each CDS in a transcript.

   gene_stats
        Table with nucleotide composition of each gene aggregated over
        all transcripts

   transcript_stats
        Table with nucleotide composition of each transcript.

   transcript_stats
        Table with amino acid composition of each protein product.

Section: geneset
----------------

Annotations derived from the ENSEMBL gene set. Annotations in
this section have been computed from the ENSEMBL gene set.
Results are in the directory :file:`geneset.dir`.

One group of :term:`bed` files outputs regions spanning whole
transcripts, genes, transcription start sites or transcription
termination sites for each of the gene sets build in
the ensembl section. These files are called
``<geneset>_<subset>_<region>.bed.gz``.

geneset
   coding, noncoding, lincrna

subset
   transcript
      regions on a per-transcript level, multiple entries
      per gene, one for each transcript
   gene
      regions on a per-gene level, one entry for each gene

regions
   region
      the complete region spanning a transcript or gene
   tss
      the transcription start site a transcript. For genes,
      it is the most upstream TSS within a gene that is reported.
   tts
      the transcription termination site of a transcript. For genes,
      it is the most downstream TTS within a gene that is reported.
   tssregion
      the region spanning all TSS within a gene

The pipeline will also compute intergenic regions for each of these
datasets called ``<geneset>_intergenic.bed.gz``.

Note that the pipeline will not compute upstream or downstream flanks
or define promotor regions as the extend of these are usually project
specific. However, these files can be easily created using bed-tools
commands taking as input the files above.

Other files in this section are:

pseudogenes.gtf.gz
   A :term:`gtf` formatted file with pseudogenes. Pseudogenes are
   either taken from the ENSEMBL annotation or processed
   transcripts with similarity to protein coding sequence. As some
   protein coding genes contain processed transcripts without an
   ORF, Pseudogenes might overlap with protein coding transcripts
   This set is not guaranteed to be complete.

numts.gtf.gz
   set of potential numts. This set is not guaranteed to be complete.

Section: gwas
-------------

Data derived from GWAS databases. The files in this section represent
regions around SNPs that have been associated with certain traits or
diseases in GWAS experiments.

gwas_catalog.bed.gz
   :term:`bed` formatted file with intervals associated with various
   traits from the `gwas catalog`_. Regions are centered around the
   listed SNPs and extended by a certain amount.

gwas_distild.bed.gz
   :term:`bed` formatted file with LD blocks associated with various traits
   from the DistilD_ database


Note that the GWAS section is only available for human.

Section: ontologies
--------------------

Data in this section are ontology assignments for genes in the ENSEMBL
geneset.

go_ensembl.tsv.gz
   table with GO assignments for genes. GO assignments are downloaded
   from ENSEMBL.

goslim_ensembl.tsv.gz
   table with GOSlim assignments for genes

go_geneontology.tsv.gz
    table with terms from geneontology.org

go_geneontology_imputed.tsv.gz
    table with terms from geneontology.org, ancestral terms imputed.

kegg
    table with imported KEGG annnotations through biomart. Note
    that KEGG through this source might be out-of-date.

Tables:
   go_ensembl_assignments
      Table with GO assignments for each gene in ENSEMBL.

   goslim_ensembl_assignments
      Table with GOSlim assignments for each gene in ENSEMBL.

   kegg_assignments
      KEGG assignments

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

For each :term:`bed`, :term:`gff` or :term:`gtf` file there is a
summary in the database called <file>_<format>_summary.  The summary
contains the number of intervals, nucleotides covered, etc. for that
particular file.

For :term:`gtf` files there is also a file with summary statistics
called <file>_gtf_stats.


Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz.

To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz
   tar -xvzf pipeline_annotations.tgz
   cd pipeline_annotations.dir
   python <srcdir>/pipeline_annotations.py make full

Code
====
"""
import sys
import shutil
import itertools
import csv
import re
import os
import glob
import collections
import pandas as pd
from ruffus import follows, transform, merge, mkdir, files, jobs_limit,\
    suffix, regex, add_inputs

import pyBigWig
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Database as Database
import CGAT.Biomart as Biomart
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineGO as PipelineGO
import CGATPipelines.PipelineUCSC as PipelineUCSC
import CGATPipelines.PipelineKEGG as PipelineKEGG
import CGAT.Intervals as Intervals


###################################################
# Pipeline configuration
###################################################
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add automatically created files to the interface.  This is required
# when the pipeline is peek'ed.  The statement below will
# add the following to the dictionary:
#
# "geneset.dir/lincrna_gene_tss.bed.gz" maps to
# "interface_geneset_lincrna_gene_tss_bed"
PARAMS.update(dict([
    ("interface_geneset_%s" %
     re.sub("[.]", "_", os.path.basename(P.snip(x, ".gz"))), x)
    for x in glob.glob('geneset.dir/*.bed.gz')]))

# Set parameter dictionary in auxilliary modules
PipelineGeneset.PARAMS = PARAMS
PipelineGO.PARAMS = PARAMS
PipelineUCSC.PARAMS = PARAMS


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    return dbh


def connectToUCSC():
    return PipelineUCSC.connectToUCSC(
        host=PARAMS["ucsc_host"],
        user=PARAMS["ucsc_user"],
        database=PARAMS["ucsc_database"])


############################################################
# Assembly
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


@follows(mkdir('assembly.dir'))
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


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS["interface_genome_tsv"])
def buildGenomeInformation(infile, outfile):
    '''
    Compute genome composition information, such as length
    and CpG density.  Uses the CGAT script `fasta2table`.

    Parameters
    ----------
    infile: str
      infile is constructed from ``PARAMS`` variable to retrieve
      the ``genome`` :term:`fasta` file

    Returns
    -------
    outfile: str
      a text file table of contigs, length and CpG density.
      The output files is GZIP compressed
    '''

    job_memory = "10G"

    statement = '''
    cat %(infile)s
    | cgat fasta2table
        --section=length
        --section=cpg
    | gzip
    > %(outfile)s
    '''
    P.run()


@P.add_doc(P.load)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(buildGenomeInformation, suffix(".tsv.gz"), ".load")
def loadGenomeInformation(infile, outfile):
    '''load genome information.'''
    P.load(infile, outfile)

##################################################################
##################################################################
##################################################################
# build G+C segmentation
##################################################################


@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"]) + ".fasta",
       PARAMS["interface_gc_segmentation_bed"])
def buildGenomeGCSegmentation(infile, outfile):
    '''
    Segments the genome into isochores - windows according to G+C
    content.  Uses `CGAT` script `fasta2bed` to generate fixed-width
    windows with their G+C content as a score.  This is then used
    as the input for `bed2bed` which merges together adjacent or
    overlapping intervals with the same number of bases into bins
    based on their score; in this case G+C content.

    Parameters
    ----------
    infile: str
      infile is constructed from `PARAMS` variable to retrieve
      the ``genome`` :term:`fasta` file

    segmentation_window_size: int
      `PARAMS` - window size to segment the genome into

    segmentation_num_bins: str
      `PARAMS` - the number of score bins to create for interval merging

    segmentation_methods: str
      `PARAMS` - method to use for merging intervals. See `bed2bed`
      documentation for details.

    Returns
    -------
    outfile: str
      :term:`BED` format file containing genome segments with similar G+C
      content.  Output file format is `BGZIP` compressed.
    '''

    statement = '''
    cgat fasta2bed
        --method=fixed-width-windows-gc
        --window-size=%(segmentation_window_size)i
        --log=%(outfile)s.log
    < %(infile)s
    | cgat bed2bed
        --method=bins
        --num-bins=%(segmentation_num_bins)s
        --binning-method=%(segmentation_method)s
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s'''

    P.run()


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"],
                    PARAMS["genome"]) + ".fasta",
       "assembly.dir/gcprofile.bed.gz")
def runGenomeGCProfile(infile, outfile):
    '''
    Uses a HMM to profile the genome based on G+C content
    and segment into windows (isochores).

    Parameters
    ----------
    infile: str
      infile is constructed from ``PARAMS`` variable to retrieve
      the ``genome`` :term:`fasta` file

    segmentation_min_length: int
      `PARAMS` - minimum length for GCProfile segmentation

    segmentation_halting_parameter: int
      `PARAMS` - GCProfile halting parameter

    Returns
    -------
    outfile: str
      :term:`BED` format file containing predicted genome
      isochores.  Outfile is `BGZIP` compressed.
    '''

    # on some cgat109 I got libstc++ error:
    # error while loading shared libraries: libstdc++.so.5
    # cannot open shared object file: No such file or directory
    to_cluster = False

    statement = '''
    cat %(infile)s
    | cgat fasta2bed
        --verbose=2
        --method=GCProfile
        --gcprofile-min-length=%(segmentation_min_length)i
        --gcprofile-halting-parameter=%(segmentation_halting_parameter)i
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s
    '''
    P.run()


@merge(runGenomeGCProfile, PARAMS["interface_gc_profile_bed"])
def buildGenomeGCProfile(infile, outfile):
    '''
    Aggregate predicted isochores with similar G+C content into
    bins.

    Parameters
    ----------
    infile: str
      output from GCProfile in :term:`BED` format. Infile is
      `BGZIP` compressed

    segmentation_num_bins: str
      `PARAMS` - the number of score bins to create for interval merging

    segmentation_methods: str
      `PARAMS` - method to use for merging intervals. See `bed2bed`
      documentation for details.

    Returns
    -------
    outfile: str
      :term:`BED` format file containing genome segments with similar G+C
      content.  Output file format is `BGZIP` compressed.
    '''
    statement = '''
    zcat %(infile)s
    | cgat bed2bed
        --method=bins
        --num-bins=%(segmentation_num_bins)s
        --binning-method=%(segmentation_method)s
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s'''

    P.run()


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

    job_memory = "10G"

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


# -----------------------------------------------------------------
# ENSEMBL gene set
@follows(mkdir('ensembl.dir'))
@files((PARAMS["ensembl_filename_gtf"], PARAMS["general_assembly_report"]), PARAMS['interface_geneset_all_gtf'])
def buildGeneSet(infiles, outfile):
    '''output sanitized ENSEMBL geneset.

    This method outputs an ENSEMBL gene set after some sanitizing steps:

    1. Chromosome names are changed to the UCSC convention.
    2. Transcripts that are not part of the chosen genome assembly
       are removed.
    3. Chromosomes that match the regular expression specified in
       the configuration file are removed.

    Arguments
    ---------
    infiles : tuple
       ENSEMBL geneset in :term:`gtf` format.
       NCBI Assembly report in `txt` format.
    outfile : string
       geneset in :term:`gtf` format.

    '''
    gtf_file, assembly_report = infiles

    statement = ['''zcat %(gtf_file)s
    | grep 'transcript_id'
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=ucsc
    --skip-missing
    --assembly-report=%(assembly_report)s
    ''']

    if PARAMS["ensembl_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --contig-pattern="%(ensembl_remove_contigs)s" ''')

    statement.append(
        '''
        | cgat gtf2gtf
        --method=set-gene_biotype-to-source
        --log=%(outfile)s.log
        | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()


@P.add_doc(PipelineGeneset.buildFlatGeneSet)
@files(buildGeneSet, PARAMS['interface_geneset_flat_gtf'])
def buildFlatGeneSet(infile, outfile):
    PipelineGeneset.buildFlatGeneSet(infile, outfile)


@P.add_doc(PipelineGeneset.loadGeneInformation)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"], "ensembl.dir/gene_info.load")
def loadGeneInformation(infile, outfile):
    '''load the transcript set.'''
    PipelineGeneset.loadGeneInformation(infile, outfile)


@P.add_doc(PipelineGeneset.loadGeneStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(mkdir('ensembl.dir'))
@files(buildFlatGeneSet, "ensembl.dir/gene_stats.load")
def loadGeneStats(infile, outfile):
    PipelineGeneset.loadGeneStats(infile, outfile)


@P.add_doc(PipelineGeneset.buildCDS)
@files(buildGeneSet,
       PARAMS["interface_geneset_cds_gtf"])
def buildCDSTranscripts(infile, outfile):
    PipelineGeneset.buildCDS(infile, outfile)


@P.add_doc(PipelineGeneset.buildExons)
@files(buildGeneSet,
       PARAMS["interface_geneset_exons_gtf"])
def buildExonTranscripts(infile, outfile):
    PipelineGeneset.buildExons(infile, outfile)


@P.add_doc(PipelineGeneset.buildCodingExons)
@files(buildGeneSet,
       PARAMS["interface_geneset_coding_exons_gtf"])
def buildCodingExonTranscripts(infile, outfile):
    PipelineGeneset.buildCodingExons(infile, outfile)


@P.add_doc(PipelineGeneset.buildNonCodingExons)
@files(buildGeneSet,
       PARAMS["interface_geneset_noncoding_exons_gtf"])
def buildNonCodingExonTranscripts(infile, outfile):
    PipelineGeneset.buildNonCodingExons(infile, outfile)


@P.add_doc(PipelineGeneset.buildLincRNAExons)
@files(buildGeneSet,
       PARAMS["interface_geneset_lincrna_exons_gtf"])
def buildLincRNAExonTranscripts(infile, outfile):
    # Jethro - some ensembl annotations contain no lincRNAs
    try:
        PipelineGeneset.buildLincRNAExons(infile, outfile)
    except Exception:
        if os.path.exists(outfile):
            assert len(IOTools.openFile(outfile).readlines()) == 0
        else:
            raise Exception("Failed to create %s" % outfile)


@P.add_doc(PipelineGeneset.loadTranscripts)
@transform((buildGeneSet,
            buildCDSTranscripts,
            buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           suffix(".gtf.gz"), "_gtf.load")
def loadTranscripts(infile, outfile):
    PipelineGeneset.loadTranscripts(infile, outfile)


@transform(buildGeneSet,
           suffix(".gtf.gz"),
           "_gtf_genome_coordinates.load")
def loadGeneCoordinates(infile, outfile):
    '''load the coordinates for each gene'''
    PipelineGeneset.loadGeneCoordinates(infile, outfile)


@P.add_doc(PipelineGeneset.loadTranscriptStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(
    ((buildExonTranscripts, "ensembl.dir/transcript_stats.load"),
     (buildCDSTranscripts, "ensembl.dir/cds_stats.load")))
def loadTranscriptStats(infile, outfile):
    PipelineGeneset.loadTranscriptStats(infile, outfile)


@jobs_limit(PARAMS.get("jobs_limit_R", 1), "R")
@follows(mkdir('ensembl.dir'))
@files(buildGeneSet, "ensembl.dir/transcript_info.load")
def downloadTranscriptInformation(infile, outfile):
    '''download information on transcripts from biomart and upload
    into database.

    This method downloads information on transcripts from the
    :term:`biomart` database and uploads it into the pipelines
    database. The columns in the mart are mapped to the following
    columns:

    * ensembl_gene_id: gene_id
    * ensembl_transcript_id: transcript_id
    * ensembl_peptide_id: protein_id
    * gene_biotype: gene_biotype
    * transcript_biotype: transcript_biotype
    * source: source
    * status: gene_status
    * transcript_status: transcript_status
    * external_gene_id: gene_name

    Only transcripts within the mart and within the supplied
    gene set are uploaded.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename with logging information. The table name
       is derived from outfile.
    ensembl_biomart_mart : PARAMS
       Biomart mart to use.
    ensembl_biomart_dataset : PARAMS
       Biomart dataset to use.
    ensembl_biomart_host : PARAMS
       Biomart host to use.
    genome : PARAMS
       Genome assembly to use. Used add missing columns
       in mart to output table.
    '''

    tablename = P.toTable(outfile)

    # use the GTF parsing approach to load the transcript information table
    PipelineGeneset.loadEnsemblTranscriptInformation(ensembl_gtf=PARAMS['ensembl_filename_gtf'],
                                                     geneset_gtf=infile,
                                                     outfile=outfile,
                                                     csvdb=PARAMS['database_name'],
                                                     set_biotype=False,
                                                     set_transcript_support=False)

    # validate: 1:1 mapping between gene_ids and gene_names
    dbh = connect()

    data = Database.executewait(dbh, """
    SELECT gene_name, count(distinct gene_id) from %(tablename)s
    GROUP BY gene_name
    HAVING count(distinct gene_id) > 1""" % locals())

    l = data.fetchall()
    if len(l) > 0:
        E.warn("there are %i gene_names mapped to different gene_ids" % len(l))
    for gene_name, counts in l:
        E.info("ambiguous mapping: %s->%i" % (gene_name, counts))

    # adding final column back into transcript_info for Drosohila and yeast
    if PARAMS["genome"].startswith("dm") or PARAMS["genome"].startswith("sac"):
        Database.executewait(
            dbh,
            '''ALTER TABLE %(tablename)s ADD COLUMN uniprot_name NULL''' %
            locals())

    P.touch(outfile)


@jobs_limit(PARAMS.get("jobs_limit_R", 1), "R")
@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"],
       "ensembl.dir/ensembl_to_entrez.load")
def downloadEntrezToEnsembl(infile, outfile):
    '''download entrez gene identifiers from biomart and upload into
    database.

    This method downloads entrez transcript identifiers from the
    :term:`biomart` database and uploads it into the pipelines
    database. The columns in the mart are mapped to the following
    columns:

    * ensembl_gene_id: gene_id
    * entrezgene: entrez_id

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename with logging information. The table name
       is derived from outfile.
    ensembl_biomart_mart : PARAMS
       Biomart mart to use.
    ensembl_biomart_dataset : PARAMS
       Biomart dataset to use.
    ensembl_biomart_host : PARAMS
       Biomart host to use.
    biomart_ensemble_gene_id : PARAMS
        Biomart attribute containing ensembl gene id
    biomart_entrez_gene_id : PARAMS
        Biomart attribute containing entrez gene id

    '''

    # SCRUM note - paramterised features being selected from biomaRt
    # in the ini file

    if not PARAMS["ensembl_biomart_mart"]:
        # skip
        P.touch(outfile)
        return None

    tablename = P.toTable(outfile)

    columns = {
        PARAMS["biomart_ensembl_gene_id"]: "gene_id",
        PARAMS["biomart_entrez_gene_id"]: "entrez_id"
        }

    data = Biomart.biomart_iterator(
        columns.keys(),
        biomart=PARAMS["ensembl_biomart_mart"],
        dataset=PARAMS["ensembl_biomart_dataset"],
        host=PARAMS["ensembl_biomart_host"])

    P.importFromIterator(
        outfile,
        tablename,
        data,
        columns=columns,
        indices=("gene_id", "entrez_id"))


@jobs_limit(PARAMS.get("jobs_limit_R", 1), "R")
@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"],
       "ensembl.dir/transcript_synonyms.load")
def downloadTranscriptSynonyms(infile, outfile):
    """download transcript synonyms from biomart and upload into database.

    This method downloads entrez transcript identifiers from the
    :term:`biomart` database and uploads it into the pipelines
    database. The columns in the mart are mapped to the following
    columns:

    * ensembl_transcript_id: transcript_id
    * external_transcript_id: transcript_name
    * refseq_mrna: refseq_id

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename with logging information. The table name
       is derived from outfile.
    ensembl_biomart_mart : PARAMS
       Biomart mart to use.
    ensembl_biomart_dataset : PARAMS
       Biomart dataset to use.
    ensembl_biomart_host : PARAMS
       Biomart host to use.
    biomart_ensemble_transcript_id : PARAMS
        Biomart attribute containing ensembl transcript id
    biomart_transcript_name : PARAMS
        Biomart attribute containing transcript name
    biomart_refseq_id : PARAMS
        Biomart attribute containing refseq ids
    """

    # SCRUM note - paramterised features being selected from biomaRt
    # in the ini file

    if not PARAMS["ensembl_biomart_mart"]:
        # skip
        P.touch(outfile)
        return None

    tablename = P.toTable(outfile)

    columns = {
        PARAMS["biomart_ensembl_transcript_id"]: "transcript_id",
        PARAMS["biomart_transcript_name"]: "transcript_name",
        PARAMS["biomart_refseq_id"]: "refseq_id"
        }

    data = Biomart.biomart_iterator(
        columns.keys(),
        biomart=PARAMS[
            "ensembl_biomart_mart"],
        dataset=PARAMS[
            "ensembl_biomart_dataset"],
        host=PARAMS["ensembl_biomart_host"])

    P.importFromIterator(
        outfile,
        tablename,
        data,
        columns=columns,
        indices=(
            "transcript_id", "transcript_name", "refseq_id"))


@P.add_doc(PipelineGeneset.buildPeptideFasta)
@follows(mkdir('ensembl.dir'))
@files(((PARAMS["ensembl_filename_pep"],
         PARAMS["interface_peptides_fasta"]), ))
def buildPeptideFasta(infile, outfile):
    PipelineGeneset.buildPeptideFasta(infile, outfile)


@P.add_doc(PipelineGeneset.buildCDNAFasta)
@follows(mkdir('ensembl.dir'))
@files(((PARAMS["ensembl_filename_cdna"],
         PARAMS["interface_cdna_fasta"]), ))
def buildCDNAFasta(infile, outfile):
    PipelineGeneset.buildCDNAFasta(infile, outfile)


@P.add_doc(PipelineGeneset.buildCDSFasta)
@follows(mkdir('ensembl.dir'))
@files((buildCDSTranscripts,
        buildPeptideFasta,),
       PARAMS["interface_cds_fasta"])
def buildCDSFasta(infiles, outfile):
    PipelineGeneset.buildCDSFasta(infiles, outfile)


@P.add_doc(PipelineGeneset.loadProteinStats)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_pep"],
       "ensembl.dir/protein_stats.load")
def loadProteinStats(infile, outfile):
    '''load the transcript set.'''
    PipelineGeneset.loadProteinStats(infile, outfile)


@merge((loadProteinStats, downloadTranscriptInformation),
       "ensembl.dir/seleno.list")
def buildSelenoList(infile, outfile):
    """export a table of seleno cysteine transcripts.

    Selenocysteine containing transcripts are identified by checking
    if their protein sequence contains ``U``.

    The table contains a single column ``transcript_id`` with ENSEMBL
    transcript identifiers as values.

    Arguments
    ---------
    infiles : list
       Unused.
    outfile : string
       Output filename in :term:`tsv` format.

    """
    # Not sure when this list is relevent or in what case it would be used - please add to documentation

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''
    SELECT DISTINCT transcript_id
    FROM transcript_info as t,
         protein_stats as p
    WHERE p.protein_id = t.protein_id AND
         p.nU > 0
    '''
    outf = open(outfile, "w")
    outf.write("transcript_id\n")
    outf.write("\n".join(
        [x[0]
         for x in Database.executewait(dbh, statement)]) + "\n")
    outf.close()


# ---------------------------------------------------------------
# geneset derived annotations
@follows(mkdir('geneset.dir'))
@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_region.bed.gz')
def buildTranscriptRegions(infile, outfile):
    """export a table of seleno cysteine transcripts.

    Selenocysteine containing transcripts are identified by checking
    if their protein sequence contains ``U``.

    The table contains a single column ``transcript_id`` with ENSEMBL
    transcript identifiers as values.

    Arguments
    ---------
    infiles : list
       Unused.
    outfile : string
       Output filename in :term:`tsv` format.

    """
    # THIS DOCUMENTATION IS NOT CORRECT - THIS NEEDS TO BE UPDATED
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_tss.bed.gz')
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_tts.bed.gz')
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
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


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
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

    statement = '''zcat %(infile)s
    | sort -k1,1 -k2,2n
    | complementBed -i stdin -g %(contigs)s
    | gzip
    > %(outfile)s'''
    P.run()


# ---------------------------------------------------------------
# UCSC derived annotations
@P.add_doc(PipelineUCSC.getRepeatsFromUCSC)
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_rna_gff"]), ))
def importRNAAnnotationFromUCSC(infile, outfile):
    """This task downloads UCSC repetetive RNA types.
    """
    # SCRUM NOTE - Why are we access ing UCSC here
    # is this a legacy thing? Andreas? Would it be better to access biomart?
    PipelineUCSC.getRepeatsFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_rnatypes"]),
        outfile=outfile,
        remove_contigs_regex=PARAMS["ensembl_remove_contigs"])


@P.add_doc(PipelineUCSC.getRepeatsFromUCSC)
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_repeats_gff"]), ))
def importRepeatsFromUCSC(infile, outfile):
    """This task downloads UCSC repeats types as identified
    in the configuration file.
    """
    # SCRUM NOTE - Why are we access ing UCSC here
    # is this a legacy thing? Andreas? Would it be better to access biomart?
    PipelineUCSC.getRepeatsFromUCSC(
        dbhandle=connectToUCSC(),
        repclasses=P.asList(PARAMS["ucsc_repeattypes"]),
        outfile=outfile)


@P.add_doc(PipelineUCSC.getCpGIslandsFromUCSC)
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_cpgislands_bed"]), ))
def importCpGIslandsFromUCSC(infile, outfile):
    '''import cpg islands from UCSC

    The repeats are stored as a :term:`bed` formatted file.
    '''
    # SCRUM NOTE - Why are we access ing UCSC here
    # is this a legacy thing? Andreas? Would it be better to access biomart?
    PipelineUCSC.getCpGIslandsFromUCSC(
        dbhandle=connectToUCSC(),
        outfile=outfile)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(importRepeatsFromUCSC, suffix(".gff.gz"), ".gff.gz.load")
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
    # SCRUM NOTE - Why are we access ing UCSC here
    # is this a legacy thing? Andreas? Would it be better to access biomart?
    load_statement = P.build_load_statement(
        tablename="repeats",
        options="--add-index=class "
        "--header-names=contig,start,stop,class")

    statement = """zcat %(infile)s
    | cgat gff2bed --set-name=class
    | grep -v "#"
    | cut -f1,2,3,4
    | %(load_statement)s
    > %(outfile)s"""
    P.run()


@transform(loadRepeats, suffix(".gff.gz.load"), ".counts.load")
def countTotalRepeatLength(infile, outfile):
    """compute genomic coverage per repeat class and load into database.

    This method computes the bases covered by each repeat class and
    uploads it into the database.
    """
    dbhandle = sqlite3.connect(PARAMS["database_name"])
    cc = dbhandle.cursor()
    statement = """DROP TABLE IF EXISTS repeat_length"""
    Database.executewait(dbhandle, statement)

    statement = """create table repeat_length as
    SELECT sum(stop-start) as total_repeat_length from repeats"""
    Database.executewait(dbhandle, statement)

    P.touch(outfile)


@P.add_doc(PipelineUCSC.getRepeatsFromUCSC)
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_allrepeats_gff"]), ))
def importAllRepeatsFromUCSC(infile, outfile):
    """This task downloads all UCSC repeats types."""
    PipelineUCSC.getRepeatsFromUCSC(dbhandle=connectToUCSC(),
                                    repclasses=None,
                                    outfile=outfile)


@follows(mkdir('ucsc.dir'))
@transform(os.path.join(PARAMS["ucsc_dir"],
                        "gbdb",
                        PARAMS["ucsc_database"],
                        "bbi",
                        "*rgMapability*.bw"),
           regex(".*rgMapabilityAlign(\d+)mer.bw"),
           add_inputs(os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"] + ".fasta")),
           r"ucsc.dir/mapability_\1.bed.gz")
def buildMapableRegions(infiles, outfile):
    '''build :term:`bed` file with mapable regions.

    Convert :term:`bigwig` data with mapability information per
    genomic position to a :term:`bed`-formatted file that lists the
    mapable regions of the genome.

    For the purpose of these tracks, a region is defined to be
    un-mapable if its maximum mapability score is less than
    0.5. Unmapable positions that are less than half the kmer size
    away from the next mapable position are designated as mapable.

    This method assumes that files use the ``CRG Alignability
    tracks``.

    UCSC says:

      The CRG Alignability tracks display how uniquely k-mer sequences
      align to a region of the genome. To generate the data, the
      GEM-mappability program has been employed. The method is
      equivalent to mapping sliding windows of k-mers (where k has been
      set to 36, 40, 50, 75 or 100 nts to produce these tracks) back to
      the genome using the GEM mapper aligner (up to 2 mismatches were
      allowed in this case). For each window, a mapability score was
      computed (S = 1/(number of matches found in the genome): S=1 means
      one match in the genome, S=0.5 is two matches in the genome, and
      so on). The CRG Alignability tracks were generated independently
      of the ENCODE project, in the framework of the GEM (GEnome
      Multitool) project.

    Arguments
    ---------
    infiles : list
       Filenames in :term:`bigwig` format with mapable data.
    outfile : string
       Output filename in :term:`bed` format with mapable regions.

    '''

    infile, fastafile = infiles
    fasta = IndexedFasta.IndexedFasta(P.snip(fastafile, ".fasta"))
    contigs = fasta.getContigSizes(with_synonyms=False)

    kmersize = int(re.search(".*Align(\d+)mer.bw", infile).groups()[0])

    E.info("creating mapable regions bed files for kmer size of %i" % kmersize)

    max_distance = kmersize // 2

    bw = pyBigWig.open(infile)

    def _iter_mapable_regions(bw, contig, size):

        min_score = PARAMS["ucsc_min_mappability"]

        # there is no iterator access, results are returned as list
        # thus proceed window-wise in 10Mb windows
        window_size = 10000000
        last_start, start = None, None

        for window_start in range(0, size, window_size):
            values = bw.intervals(contig, window_start, window_start + window_size)
            if values is None:
                continue

            for this_start, this_end, value in values:
                if value < min_score:
                    if start:
                        yield start, this_start
                    start = None
                else:
                    if start is None:
                        start = this_start

        if start is not None:
            yield start, this_end

    outf = IOTools.openFile(outfile, "w")

    for contig, size in contigs.items():

        last_start, last_end = None, None
        for start, end in _iter_mapable_regions(bw, contig, size):
            if last_start is None:
                last_start, last_end = start, end
            if start - last_end >= max_distance:
                outf.write("%s\t%i\t%i\n" % (contig, last_start, last_end))
                last_start = start

            last_end = end

        if last_start is not None:
            outf.write("%s\t%i\t%i\n" % (contig, last_start, last_end))

    outf.close()


@transform(buildMapableRegions, suffix(".bed.gz"),
           ".filtered.bed.gz")
def filterMapableRegions(infile, outfile):
    """remove small windows from a mapability track.

    Too many fragmented regions will cause gat to fail as it fragments
    the workspace in a GAT analysis into too many individual segments.

    The filtering works by merging all segments that are within
    mapability_merge_distance and removing all those that are larger
    than mapabpility_min_segment_size.

    Arguments
    ---------
    infile : string
       Input filename in :term:`bed` format.
    outfile : string
       Output filename in :term:`bed` format with mapable regions.
    mapability_merge_distance : int
       see :term:`PARAMS`
    mapability_min_segment_size : int
       see :term:`PARAMS`

    """

    statement = '''
    mergeBed -i %(infile)s -d %(mapability_merge_distance)i
    | awk '$3 - $2 >= %(mapability_min_segment_size)i'
    | gzip
    > %(outfile)s
    '''

    P.run()


# ---------------------------------------------------------------
# GWAS data
if PARAMS["genome"].startswith("hg"):

    @follows(mkdir('gwas.dir'))
    @merge(None, "gwas.dir/gwascatalog.txt")
    def downloadGWASCatalog(infile, outfile):
        '''
        Download the GWAS catalog data for the human genome

        Parameters
        ----------
        infile: None
          an unused variable required by Ruffus

        Returns
        -------
        outfile: str
          an `excel` file containing the human genome GWAS catalog
        '''

        if os.path.exists(outfile):
            os.remove(outfile)
            # MM: this is hard-coded - the URL can (and has) changed, so
            # this should be defined in the pipeline config file
            # AH: Moved to EBI, download needs to be updated
        statement = '''curl https://www.genome.gov/admin/gwascatalog.txt
        | sed 's/[\d128-\d255]//g'
        > %(outfile)s'''
        P.run()

    @merge(downloadGWASCatalog, PARAMS["interface_gwas_catalog_bed"])
    def buildGWASCatalogTracks(infile, outfile):
        '''
        Convert the GWAS catalog entries to :term:`BED` format.

        Parameters
        ----------
        infile: str
          an `excel` format file of GWAS catalog entries

        genome_dir: str
          PARAMS - directory containing the indexed :term:`FASTA` genome
          files

        genome: str
          PARAMS - indexed genome build to use

        gwas_extension: int
          PARAMS - size in bp to extend region around each GWAS catalog entry

        Returns
        -------
        outfile: str
          :term:`BED` format file of GWAS catalog entries
        '''

        reader = csv.DictReader(IOTools.openFile(infile),
                                dialect="excel-tab")

        tracks = collections.defaultdict(lambda: collections.defaultdict(list))

        fasta = IndexedFasta.IndexedFasta(
            os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"))
        contigsizes = fasta.getContigSizes()
        c = E.Counter()

        for row in reader:
            c.input += 1
            contig, pos, snp, disease = row['Chr_id'], row[
                'Chr_pos'], row['SNPs'], row['Disease/Trait']

            # skip SNPs on undefined contigs
            if contig not in contigsizes:
                c.no_contig += 1
                continue

            if snp == "NR":
                c.skipped += 1
                continue

            if pos == "":
                c.no_pos += 1
                continue

            # translate chr23 to X
            if contig == "23":
                contig = "X"

            contig = "chr%s" % contig

            try:
                tracks[disease][contig].append(int(pos))
            except ValueError:
                print(row)
            c.output += 1

        E.info(c)

        extension = PARAMS["gwas_extension"]

        c = E.Counter()
        outf = IOTools.openFile(outfile, "w")
        for disease, pp in tracks.items():

            for contig, positions in pp.items():
                contigsize = contigsizes[contig]
                regions = [(max(0, x - extension),
                            min(contigsize, x + extension))
                           for x in positions]

                regions = Intervals.combine(regions)
                c[disease] += len(regions)

                for start, end in regions:
                    outf.write("%s\t%i\t%i\t%s\n" %
                               (contig, start, end, disease))

        outf.close()

        outf = IOTools.openFile(outfile + ".log", "w")
        outf.write("category\tcounts\n%s\n" % c.asTable())
        outf.close()

    @follows(mkdir('gwas.dir'))
    @merge(None, "gwas.dir/gwas_distild.log")
    def downloadDistiLD(infile, outfile):
        '''
        Download GWAS data from the DistiLD database.

        Parameters
        ----------
        infile: None
          an unused variable required by Ruffus

        Returns
        -------
        outfile: str
          two text files are output that contain SNP LD blocks with
          gene annotations and SNP IDs, and SNP IDs with GWAS
          associations and linked ICD10 codes
        '''

        track = P.snip(outfile, ".log")
        of = track + "_snps.tsv.gz"
        if os.path.exists(of):
            os.remove(of)
        statement = \
            '''wget http://distild.jensenlab.org/snps.tsv.gz
            -O %(of)s'''
        P.run()

        of = track + "_lds.tsv.gz"
        if os.path.exists(of):
            os.remove(of)
        statement = \
            '''wget http://distild.jensenlab.org/lds.tsv.gz
            -O %(of)s'''
        P.run()

        P.touch(outfile)

    @merge(downloadDistiLD, PARAMS["interface_gwas_distild_bed"])
    def buildDistiLDTracks(infile, outfile):
        '''
        Build :term:`BED` tracks from entries in the DistiLD database
        of disease/trait associations

        Parameters
        ----------
        infile: str
          the log file from the downloading DistiLD database files

        genome_dir: str
          PARAMS - directory containing the indexed :term:`FASTA` genome
          files

        genome: str
          PARAMS - indexed genome build to use

        Returns
        -------
        outfile: str
          :term:`BED` format file containing disease associated SNPs
          and their associated trait(s)
        '''

        track = P.snip(infile, ".log")
        intervals = []

        fasta = IndexedFasta.IndexedFasta(
            os.path.join(PARAMS["genome_dir"],
                         PARAMS["genome"] + ".fasta"))
        contigsizes = fasta.getContigSizes()

        c = E.Counter()
        for line in IOTools.openFile(track + "_snps.tsv.gz"):
            pubmed_id, rs, pvalue, block, ensgenes, short, icd10 = line[
                :-1].split("\t")
            c.input += 1
            try:
                contig, start, end = re.match(
                    "(\S+):(\d+)-(\d+)", block).groups()
            except AttributeError:
                E.warn("parsing error for %s" % block)
                c.errors += 1
                continue

            # skip SNPs on undefined contigs
            if contig not in contigsizes:
                c.no_contig += 1
                continue

            intervals.append((contig, int(start), int(end), short))
            c.parsed += 1

        intervals.sort()
        outf = IOTools.openFile(outfile, "w")
        cc = E.Counter()
        for k, x in itertools.groupby(intervals, key=lambda x: x):
            outf.write("%s\t%i\t%i\t%s\n" % k)
            c.output += 1
            cc[k[3]] += 1
        outf.close()
        E.info(c)

        outf = IOTools.openFile(outfile + ".log", "w")
        outf.write("category\tcounts\n%s\n" % cc.asTable())
        outf.close()

    @follows(buildGWASCatalogTracks, buildDistiLDTracks)
    def _gwas():
        pass
else:
    @files(((None, None),))
    def _gwas(infile, outfile):
        pass


# ---------------------------------------------------------------
# Ontologies
# SCRUM NOTES - ARE THESE LEGACY FROM OLD ENRICHMENT PIPELINE?
# Can we remove them to streamline the pipeline - its failing here
# on KEGG and on GO - No don't remove just make it work

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
    PipelineGO.createGOFromENSEMBL(infile, outfile)


@P.add_doc(PipelineGO.createGOSlimFromENSEMBL)
@transform(createGO,
           regex("(.*)"),
           PARAMS["interface_goslim_ensembl"])
def createGOSlim(infile, outfile):
    '''
    Downloads GO slim annotations from ensembl
    '''
    PipelineGO.createGOSlimFromENSEMBL(infile, outfile)


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


@P.add_doc(PipelineGO.buildGOPaths)
@transform(createGO, suffix(".tsv.gz"), ".paths")
def buildGOPaths(infile, outfile):
    '''compute a file with paths of each GO term to the ancestral node.'''
    infile = P.snip(infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOPaths(infile, outfile)


@P.add_doc(PipelineGO.buildGOTable)
@transform(createGO, suffix(".tsv.gz"), ".desc.tsv")
def buildGOTable(infile, outfile):
    '''build a simple table with GO descriptions in obo.'''
    infile = P.snip(infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOTable(infile, outfile)


@P.add_doc(P.load)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(buildGOTable, suffix(".tsv"), ".load")
def loadGOTable(infile, outfile):
    '''load GO descriptions into database.'''
    P.load(infile, outfile)


@P.add_doc(PipelineGO.createGOFromGeneOntology)
@follows(mkdir('ontologies.dir'),
         downloadTranscriptInformation, loadGOAssignments)
@files([(None, PARAMS["interface_go_geneontology"]), ])
def createGOFromGeneOntology(infile, outfile):
    '''build GO assignments from GeneOntology.org'''
    PipelineGO.createGOFromGeneOntology(infile, outfile)


@P.add_doc(PipelineGO.imputeGO)
@transform(createGOFromGeneOntology,
           suffix(".tsv.gz"),
           add_inputs(buildGOPaths),
           PARAMS["interface_go_geneontology_imputed"])
def imputeGO(infiles, outfile):
    '''imput ancestral GO terms for each gene based on
    derived GO terms.
    '''
    PipelineGO.imputeGO(infiles[0], infiles[1], outfile)

# THIS IS CURRRENTLY FAILYING - NEED TO CHECK R CODE
# AND FIX

# I have fixed it in a commit to cgat/CGAT/Biomart.py - KB


@jobs_limit(PARAMS.get("jobs_limit_R", 1), "R")
@P.add_doc(PipelineKEGG.importKEGGAssignments)
@follows(mkdir('ontologies.dir'))
@files(None, PARAMS['interface_kegg'])
def importKEGGAssignments(infile, outfile):
    '''
    Imports the KEGG annotations from the R KEGG.db package

    Note that since KEGG is no longer
    publically availible, this is not up-to-date and maybe removed
    from bioconductor in future releases

    Entrez IDs are downloaded from Biomart
    Corresponding KEGG IDs are downloaded from KEGG.db using
    KEGGEXTID2PATHID then translated to path names using
    KEGGPATHID2NAME.
    '''

    biomart_dataset = PARAMS["KEGG_dataset"]
    mart = PARAMS["KEGG_mart"]

    # Possibly this should use the same biomart version as the rest of the
    # pipeline by calling ensembl_biomart_host instead of
    # KEGG_host from PARAMS KB
    host = PARAMS["KEGG_host"]
    PipelineKEGG.importKEGGAssignments(outfile, mart, host, biomart_dataset)


@P.add_doc(P.load)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(importKEGGAssignments, suffix(".tsv.gz"), "_assignments.load")
def loadKEGGAssignments(infile, outfile):

    P.load(infile, outfile, options="-i gene_id -i kegg_id --allow-empty-file")


# ---------------------------------------------------------------
# Enrichment analysis
@P.add_doc(PipelineGeneset.annotateGenome)
@follows(mkdir('enrichment.dir'))
@files(buildGeneSet, PARAMS['interface_annotation_gff'])
def annotateGenome(infile, outfile):
    """This task only considers protein coding genes as
    processed_transcripts tend to cover larger genomic regions and
    often overlap between adjacent protein coding genes.

    """
    # This could case problems if source collumn have changed
    # need to add in a check that this is acessing the right info
    # maybe make more explicit so that it can know the right gtf attribute to access
    # Add in ability to output some stats on how many annotations, how many are pt-coding
    PipelineGeneset.annotateGenome(infile,
                                   outfile,
                                   only_proteincoding=True)


@P.add_doc(PipelineGeneset.annotateGeneStructure)
@follows(mkdir('enrichment.dir'))
@files(buildGeneSet, PARAMS['interface_genestructure_gff'])
def annotateGeneStructure(infile, outfile):
    """This task only considers protein coding genes as
    processed_transcripts tend to cover larger genomic regions and
    often overlap between adjacent protein coding genes.

    """
    PipelineGeneset.annotateGeneStructure(infile,
                                          outfile,
                                          only_proteincoding=True)


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


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_tssterritories_gff"])
def buildTSSTerritories(infile, outfile):
    """build TSS territories from protein coding genes.

    The tss territory of a gene is defined as a region centered aronud
    the TSS. If the territories of two genes overlap, they are
    resolved at the mid-point between the two adjacent genes.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gff` format.
    enrichment_territories_radius : int
       see :term:`PARAMS`
    """
    statement = '''
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=filter
    --filter-method=proteincoding
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=filter
    --filter-method=representative-transcript
    --log=%(outfile)s.log
    | cgat gtf2gtf --method=sort --sort-order=position
    | cgat gtf2gff
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --territory-extension=%(enrichment_territories_radius)s
    --method=tss-territories
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=filter --filter-method=longest-gene --log=%(outfile)s.log
    | gzip
    > %(outfile)s '''

    P.run()


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_greatdomains_gff"])
def buildGREATRegulatoryDomains(infile, outfile):
    """build GREAT regulatory domains.

    Each TSS in a gene is associated with a basal region. The basal
    region is then extended upstream to the basal region of the
    closest gene, but at most by a certain radius. In the case of
    overlapping genes, the extension is towards the next
    non-overlapping gene.

    This is the "basal plus extension" rule in GREAT. Commonly used
    are 5+1 with 1 Mb extension.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gff` format.
    enrichment_great_radius : int
       see :term:`PARAMS`
    enrichment_great_upstream : int
       see :term:`PARAMS`
    enrichment_great_downstream : int
       see :term:`PARAMS`

    """

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter
    --filter-method=proteincoding
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=filter --filter-method=representative-transcript
    --log=%(outfile)s.log
    | cgat gtf2gff
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --method=great-domains
    --territory-extension=%(enrichment_great_radius)s
    --upstream-extension=%(enrichment_great_upstream)i
    --downstream-extension=%(enrichment_great_downstream)i
    | gzip
    > %(outfile)s '''

    P.run()


@P.add_doc(PipelineGeneset.buildGenomicContext)
@follows(mkdir('enrichment.dir'))
@merge((importRepeatsFromUCSC,
        importRNAAnnotationFromUCSC,
        buildGeneSet,
        buildFlatGeneSet,
        importCpGIslandsFromUCSC,
        createGO,
        ),
       PARAMS["interface_genomic_context_bed"])
def buildGenomicContext(infiles, outfile):
    PipelineGeneset.buildGenomicContext(infiles, outfile)
    # Scrum notes
    # This needs some attention - check the output of this between the builds
    # are all the collumn headers the same, are there similar numbers in each one
    # does this have overlapping contexts or can a region have only one context

    # This feeds down to context stats - this also needs attention after this step has been verified
    # make sure there are stats on this table in some part of the report

    # HEre are the stats - check these are reasonable and in report


@transform(buildGenomicContext, suffix(".bed.gz"), ".tsv")
def buildGenomicContextStats(infile, outfile):
    """compute overlap between annotations in a :term:`bed` file.

    This method splits a :term:`bed` formatted file by its fourth
    column, the feature name. It then computes the individual :term:`bed`
    formatted files with :doc:`diff_bed`.

    Arguments
    ---------
    infiles : string
        Input filename of :term:`bed` formatted file with annotations.
    outfile : string
        Output filename in :term:`tsv` format.
    """

    tmpdir = P.getTempDir(".")

    statement = '''zcat %(infile)s
    | cgat split_file
        --pattern-output=%(tmpdir)s/%%s.bed
        --column=4
    > %(outfile)s.log
    '''

    P.run()

    statement = '''
    cgat diff_bed
       %(tmpdir)s/*.bed
    > %(outfile)s
    '''
    P.run()

    shutil.rmtree(tmpdir)


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
        outfiles=outfiles)


@P.add_doc(PipelineGeneset.buildPseudogenes)
@files((buildGeneSet,
        buildPeptideFasta),
       PARAMS["interface_pseudogenes_gtf"])
def buildPseudogenes(infile, outfile):
    dbh = connect()
    PipelineGeneset.buildPseudogenes(infile, outfile, dbh)


@P.add_doc(PipelineGeneset.buildNUMTs)
@follows(mkdir('geneset.dir'))
@files((None,),
       PARAMS["interface_numts_gtf"])
def buildNUMTs(infile, outfile):
    PipelineGeneset.buildNUMTs(infile, outfile)

# --------------------------------------------
# Below is a collection of functions that are
# currently inactivated.

# This is all legacy - sebastian says this not appropriate programming
# behavoir :'(
if 0:
    ############################################################
    ############################################################
    ############################################################
    # get UCSC tables
    ############################################################
    def getUCSCTracks(infile=PARAMS["filename_ucsc_encode"]):
        '''return a list of UCSC tracks from infile.'''
        tables = []
        with open(infile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                tablename = line[:-1].strip()
                if tablename == "":
                    continue
                tables.append(tablename)
        return tables

    ############################################################
    ############################################################
    ############################################################
    # import UCSC encode tracks
    ############################################################
    @posttask(touch_file("ucsc_encode.import"))
    @files(PARAMS["filename_ucsc_encode"], "ucsc_encode.import")
    def importUCSCEncodeTracks(infile, outfile):

        dbhandle = sqlite3.connect(PARAMS["database_name"])

        cc = dbhandle.cursor()
        tables = set(
            [x[0] for x in cc.executewait(
                dbhandle,
                "SELECT name FROM sqlite_master WHERE type='table'")])
        cc.close()

        for tablename in getUCSCTracks(infile):
            if tablename in tables:
                E.info("skipping %(tablename)s - already exists" % locals())
                continue

            load_statement = P.build_load_statement(tablename)

            E.info("importing %(tablename)s" % locals())

            statement = '''
            mysql --user=genome --host=genome-mysql.cse.ucsc.edu
            -A -B -e "SELECT * FROM %(tablename)s" %(ucsc_database)s
            | %(load_statement)s
            >> %(outfile)s
            '''
            P.run()

    ############################################################
    ############################################################
    ############################################################
    # export UCSC encode tracks as bed
    ############################################################
    @transform(importUCSCEncodeTracks, suffix(".import"), ".bed")
    def exportUCSCEncodeTracks(infile, outfile):

        dbhandle = sqlite3.connect(PARAMS["database_name"])

        outs = open(outfile, "w")
        for tablename in getUCSCTracks():
            outs.write("track name=%s\n" % tablename)

            cc = dbhandle.cursor()
            statement = """SELECT chrom, chrostart, chroend FROM %s
            ORDER by chrom, chrostart""" % (
                tablename)
            cc.executewait(dbhandle, statement)
            for contig, start, end in cc:
                outs.write("%s\t%i\t%i\n" % (contig, start, end))
        outs.close()


@transform("*/*.gff.gz",
           suffix(".gff.gz"),
           ".gffsummary.tsv.gz")
def buildGFFSummary(infile, outfile):
    """summarize genomic coverage of a :term:`gff` formatted file.

    Arguments
    ---------
    infile : string
        Input filename of :term:`gff` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.

    """
    statement = '''zcat %(infile)s
    | cgat gff2coverage
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.bed.gz",
           suffix(".bed.gz"),
           ".bedsummary.tsv.gz")
def buildBedSummary(infile, outfile):
    """summarize genomic coverage of a :term:`bed` formatted file.

    The coverage is computed per contig.

    Arguments
    ---------
    infile : string
        Input filename of :term:`bed` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.

    """
    statement = '''zcat %(infile)s
    | cgat bed2stats
    --aggregate-by=contig
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/genomic_context.bed.gz",
           suffix(".bed.gz"),
           ".bednamesummary.tsv.gz")
def buildBedNameSummary(infile, outfile):
    """summarize genomic coverage of a :term:`bed` formatted file.

    The coverage is computed per annotation (column 4) in the
    :term:`bed` file.

    Arguments
    ---------
    infile : string
        Input filename of :term:`bed` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.

    """
    statement = '''zcat %(infile)s
    | cgat bed2stats
    --aggregate-by=name
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gtf.gz",
           suffix(".gtf.gz"),
           ".gtfsummary.tsv.gz")
def buildGTFSummary(infile, outfile):
    """summarize genomic coverage of a :term:`gtf` formatted file.

    Arguments
    ---------
    infile : string
        Input filename of :term:`gtf` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.
    """

    statement = '''zcat %(infile)s
    | cgat gff2coverage
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gtf.gz",
           suffix(".gtf.gz"),
           ".gtfstats.tsv.gz")
def buildGTFStats(infile, outfile):
    """summarize stats of a :term:`gtf` formatted file.

    The statistics are number of genes, transcripts, etc.

    Arguments
    ---------
    infile : string
        Input filename of :term:`gtf` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.
    """
    statement = '''zcat %(infile)s
    | cgat gff2stats
    --is-gtf
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gff.gz",
           suffix(".gff.gz"),
           ".gffstats.tsv.gz")
def buildGFFStats(infile, outfile):
    """summarize stats of a :term:`gff` formatted file.

    The statistics are number of contigs, strands
    features and sources.

    Arguments
    ---------
    infile : string
        Input filename of :term:`gff` formatted file.
    outfile : string
        Output filename in :term:`tsv` format.
    """
    statement = '''zcat %(infile)s
    | cgat gff2stats
    | gzip > %(outfile)s
    '''
    P.run()


@merge(buildGTFStats, 'gtf_stats.load')
def loadGTFStats(infiles, outfile):
    """load summary data into database."""
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.*).tsv.gz",
                         options="--allow-empty")


@merge(buildGFFStats, 'gff_stats.load')
def loadGFFStats(infiles, outfile):
    """load summary data into database."""
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.*).tsv.gz",
                         options="--allow-empty")


@transform((buildGFFSummary,
            buildBedSummary,
            buildBedNameSummary,
            buildGTFSummary),
           suffix(".tsv.gz"),
           ".load")
def loadIntervalSummary(infile, outfile):
    """load summary data into database."""
    P.load(infile, outfile, options='--allow-empty-file')


##################################################################
# Primary targets
@follows(buildContigSizes,
         buildContigBed,
         buildUngappedContigBed,
         loadGenomeInformation,
         buildGenomeGCProfile,
         buildCpGBed)
def assembly():
    """convenience target : assembly derived annotations"""
    pass


@follows(buildGeneSet,
         loadTranscripts,
         downloadTranscriptInformation,
         loadGeneStats,
         loadTranscriptStats,
         loadGeneInformation,
         loadGeneCoordinates,
         downloadEntrezToEnsembl,
         downloadTranscriptSynonyms,
         buildExonTranscripts,
         buildCodingExonTranscripts,
         buildNonCodingExonTranscripts,
         buildPseudogenes,
         buildNUMTs,
         buildSelenoList,
         )
def ensembl():
    """convenience target : ENSEMBL geneset derived annotations"""
    pass


@follows(buildPeptideFasta,
         buildCDSFasta,
         buildCDNAFasta)
def fasta():
    """convenience target : sequence collections"""
    pass


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
    pass


@follows(importRepeatsFromUCSC,
         importRNAAnnotationFromUCSC,
         importCpGIslandsFromUCSC,
         loadRepeats,
         countTotalRepeatLength)
def ucsc():
    """convenience target : UCSC derived annotations"""
    pass


# annotation targets are only intrinsic data sets
# based on the genome and the gene set
@follows(buildGeneTerritories,
         buildTSSTerritories,
         buildGREATRegulatoryDomains,
         annotateGeneStructure,
         annotateGenome)
def annotations():
    """convenience target : gene based annotations"""
    pass


# enrichment targets include extrinsic data sets such
# as GO, UCSC, etc.
@follows(buildGenomicContext,
         buildGenomicContextStats,
         buildGenomicFunctionalAnnotation)
def enrichment():
    """convenience target : annotations for enrichment analysis"""
    pass


@follows(loadGOAssignments,
         loadKEGGAssignments)
def ontologies():
    """convenience target : ontology information"""
    pass


@follows(_gwas)
def gwas():
    """convenience target : import GWAS data"""
    pass


@follows(loadGTFStats,
         loadGFFStats,
         loadIntervalSummary)
def summary():
    '''convenience target : summary'''
    pass


# @follows(calculateMappability, countMappableBases,
#          loadMappableBases, splitMappabiliyFileByContig,
#          countMappableBasesPerContig, loadMappableBasesPerContig)
# def gemMappability():
#     '''Count mappable bases in genome'''
#     pass
# taken out gemMappability as not fully configured


@follows(assembly,
         ensembl,
         ucsc,
         geneset,
         # fasta,   # AH disabled for now, peptides2cds missing
         ontologies,
         annotations,
         enrichment,
         gwas)
def full():
    '''build all targets - note: run summary separately afterwards.'''
    pass

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


@follows(mkdir("report"), summary)
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"), summary)
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
