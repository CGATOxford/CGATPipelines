'''PipelineGeneset.py - Tasks for processing gene sets
======================================================

Most of this tasks take a geneset (.gtf.gz) from ENSEMBL as input.

As of ENSEMBL release 75 the gtf file contains both transcripts but
also untranscribed features such as pseudo genes, for example::

   1       pseudogene      gene    11869   14412   .       +       .       gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";

Reference
---------

'''

import os
import collections
import sqlite3

import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta

# When importing this module, set PARAMS to your parameter
# dictionary
PARAMS = {}

ENSEMBL_INFO = collections.namedtuple(
    "ENSEMBLINFO", "species gene_prefix transcript_prefix")

# Map of UCSC genome prefixes to ENSEMBL gene sets
MAP_UCSC2ENSEMBL = {
    'hg': ENSEMBL_INFO._make(('Homo_sapiens',
                              'ENSG',
                              'ENST')),
    'mm': ENSEMBL_INFO._make(('Mus_musculus',
                              'ENSMUSG',
                              'ENSMUST')),
    'rn': ENSEMBL_INFO._make(('Rattus_norvegicus',
                              'ENSRNOG',
                              'ENSRNOT')),
}


def mapUCSCToEnsembl(genome):
    '''map the name of a UCSC genome (hg19, mm10) to
    ENSEMBL URLs.'''
    prefix = genome[:2]
    return MAP_UCSC2ENSEMBL[prefix]


def annotateGenome(infile, outfile,
                   only_proteincoding=False):
    '''annotate genomic regions with reference gene set.

    The method applies the following filters to an ENSEMBL gene set:

    * Select transcribed features, i.e., those entries that contain a
      ``transcript_id``.

    * Merge overlapping exons from different transcripts within a
      gene.

    * In case of overlapping genes, take the longest gene in genomic
      coordinates.

    The resultant gene set is then converted to genomic annotations
    such as exonic, intronic, intergenic. For more information, see
    documentation for the script :mod:`gtf2gff.py` under the option
    ``--method=genome``.

    Arguments
    ---------

    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gff` format.
    only_proteincoding : bool
       If True, only consider protein coding genes.
    '''

    method = "genome"

    if only_proteincoding:
        filter_cmd = """cgat gtf2gtf
        --method=filter --filter-method=proteincoding""" % PARAMS
    else:
        filter_cmd = "cat"

    statement = """
    zcat %(infile)s
    | %(filter_cmd)s
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | cgat gtf2gtf
    --method=set-source-to-transcript_biotype
    | cgat gtf2gtf
    --method=merge-exons
    --mark-utr
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=filter --filter-method=longest-gene
        --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort --sort-order=position
    | cgat gtf2gff
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --flank-size=%(enrichment_genes_flank)s
    --method=%(method)s
    | gzip
    > %(outfile)s
    """
    P.run()


def annotateGeneStructure(infile, outfile,
                          only_proteincoding=False):
    """annotate genomic regions with gene structure.

    The method applies the following filters to an ENSEMBL gene set:

    * Select transcribed features, i.e., those entries that contain a
      ``transcript_id``.

    * If there are multiple transcripts per gene, take a
      representative transcript. See :mod:`gtf2gtf` for the definition
      of the representative transcript.

    * In case of overlapping genes, take the longest gene in genomic
      coordinates.

    The resultant gene set is then converted to genomic annotations
    such as first_exon, first_intron, .... For more information, see
    documentation for the script :mod:`gtf2gff.py` under the option
    ``--method=genes``.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gff` format.
    only_proteincoding : bool
       If True, only consider protein coding genes.

    """

    if only_proteincoding:
        filter_cmd = """cgat gtf2gtf
        --method=filter --filter-method=proteincoding""" % PARAMS
    else:
        filter_cmd = "cat"

    method = "genes"

    statement = """
    gunzip
    < %(infile)s
    | %(filter_cmd)s
    | awk '$3 == "exon"'
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | cgat gtf2gtf
    --method=filter --filter-method=representative-transcript
    | cgat gtf2gtf
    --method=filter --filter-method=longest-gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort --sort-order=position
    | cgat gtf2gff
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --flank-size=%(enrichment_genestructures_flank)i
    --flank-increment-size=%(enrichment_genestructures_increment)i
    --method=%(method)s
    --gene-detail=exons
    | gzip
    > %(outfile)s
    """
    P.run()


def buildFlatGeneSet(infile, outfile):
    '''build a flattened gene set.

    All transcripts in a gene are merged into a single transcript by
    combining overlapping exons.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''
    # sort by contig+gene, as in refseq gene sets, genes on
    # chr_random might contain the same identifier as on chr
    # and hence merging will fail.
    # --permit-duplicates is set so that these cases will be
    # assigned new merged gene ids.

    statement = """gunzip
    < %(infile)s
    | awk '$3 == "exon"'
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort
    --sort-order=contig+gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=merge-exons
    --permit-duplicates
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=set-transcript-to-gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort
    --sort-order=position+gene
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
        """
    P.run()


def buildProteinCodingGenes(infile, outfile):
    '''build a proctein coding gene set from an ENSEMBL gene set.

    The method applies the following filters to an ENSEMBL gene set:

    * Select protein coding features.

    * Remove features not on the reference genome that has been chosen.

    * Merge overlapping exons from different transcripts within a
      gene.

    * In case of overlapping genes, take the longest gene in genomic
      coordinates.

    * Keep only features called ``exon`` in the GTF file.

    * Set the ``transcript_id`` to the ``gene_id``

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''

    # sort by contig+gene, as in refseq gene sets, genes on
    # chr_random might contain the same identifier as on chr
    # and hence merging will fail.
    # --permit-duplicates is set so that these cases will be
    # assigned new merged gene ids.
    statement = """zcat %(infile)s
    | cgat gtf2gtf
        --method=filter
        --filter-method=proteincoding
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=contig+gene
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    | cgat gtf2gtf
    --method=merge-exons
    --permit-duplicates
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=filter --filter-method=longest-gene
    --log=%(outfile)s.log
    | awk '$3 == "exon"'
    | cgat gtf2gtf
    --method=set-transcript-to-gene
    --log=%(outfile)s.log
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | gzip
    > %(outfile)s
    """
    P.run()


def loadGeneInformation(infile, outfile, only_proteincoding=False):
    '''load gene-related attributes from :term:`gtf` file into database.

    This method takes transcript-associated features from an
    :term:`gtf` file and collects the gene-related attributes in the
    9th column of the gtf file, ignoring exon_id, transcript_id,
    transcript_name, protein_id and exon_number.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename, contains logging information. The
       table name is derived from the filename of outfile.
    only_proteincoding : bool
       If True, only consider protein coding genes.

    '''

    job_memory = "4G"
    table = P.toTable(outfile)

    if only_proteincoding:
        filter_cmd = """cgat gtf2gtf
        --method=filter --filter-method=proteincoding""" % PARAMS
    else:
        filter_cmd = "cat"

    load_statement = P.build_load_statement(
        table,
        options="--add-index=gene_id "
        "--add-index=gene_name"
        "--map=gene_name:str")

    statement = '''
    zcat %(infile)s
    | %(filter_cmd)s
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | cgat gtf2tsv
    --attributes-as-columns --output-only-attributes -v 0
    | python %(toolsdir)s/csv_cut.py
    --remove exon_id transcript_id transcript_name protein_id exon_number
    | %(pipeline_scriptsdir)s/hsort 1
    | uniq
    | %(load_statement)s
    > %(outfile)s'''

    P.run()


def loadTranscriptInformation(infile, outfile,
                              only_proteincoding=False):
    '''load transcript-related attributes from :term:`gtf` file into database.

    This method takes transcript-associated features from an
    :term:`gtf` file and collects the gene-related attributes in the
    9th column of the gtf file, ignoring exon_id and exon_number.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename, contains logging information. The
       table name is derived from the filename of outfile.
    only_proteincoding : bool
       If True, only consider protein coding genes.

    '''
    table = P.toTable(outfile)

    if only_proteincoding:
        filter_cmd = """cgat gtf2gtf
        --method=filter --filter-method=proteincoding""" % PARAMS
    else:
        filter_cmd = "cat"

    load_statement = P.build_load_statement(
        table,
        options="--add-index=gene_id "
        "--add-index=gene_name"
        "--add-index=protein_id"
        "--add-index=transcript_id"
        "--map=gene_name:str")

    statement = '''zcat < %(infile)s
    | awk '$3 == "CDS"'
    | grep "transcript_id"
    | cgat gtf2gtf
    --method=sort --sort-order=gene+transcript
    | cgat gtf2tsv
    --attributes-as-columns --output-only-attributes -v 0
    | python %(toolsdir)s/csv_cut.py --remove exon_id exon_number
    | %(pipeline_scriptsdir)s/hsort 1 | uniq
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


def buildCDNAFasta(infile, outfile):
    '''index an ENSEMBL cdna FASTA file

    The descriptions in the fasta file are truncated at the
    first space to contain only the sequence identifier.

    Arguments
    ---------
    infile : string
        ENSEMBL ``.cdna.fa.gz`` file in :term:`fasta` format
    outfile : string
        indexed file in :term:`fasta` format
    '''
    dbname = outfile[:-len(".fasta")]

    # perl statement to truncate ids from ENSPXXXX.1 to ENSPXXX.
    # New notation introduced around Ensembl release 80
    # Necessary for comparing to gtf work and backward compatibility.
    statement = '''gunzip
    < %(infile)s
    | perl -p -e 'if ("^>") { s/[\z\.].*//};'
    | cgat index_fasta
       --force-output
    %(dbname)s -
    > %(dbname)s.log
    '''

    P.run()


def buildPeptideFasta(infile, outfile):
    '''index an ENSEMBL peptide FASTA file

    The descriptions in the fasta file are truncated at the
    first space to contain only the sequence identifier.

    Arguments
    ---------
    infile : string
        ENSEMBL ``.pep.all.fa.gz`` file in :term:`fasta` format
    outfile : string
        indexed file in :term:`fasta` format
    '''
    dbname = outfile[:-len(".fasta")]

    # perl statement to truncate ids from ENSPXXXX.1 to ENSPXXX.
    # New notation introduced around Ensembl release 80
    # Necessary for comparing to gtf work and backward compatibility.
    statement = '''gunzip
    < %(infile)s
    | perl -p -e 'if ("^>") { s/[\z\.].*//};'
    | cgat index_fasta
       --force-output
    %(dbname)s -
    > %(dbname)s.log
    '''

    P.run()


def loadPeptideSequences(infile, outfile):
    '''load ENSEMBL peptide file into database

    This method removes empty sequences (see for example
    transcript:ENSMUST00000151316, ENSMUSP00000118372)

    The created table contains the columns ``protein_id``, ``length``
    and ``sequence``.

    Arguments
    ---------
    infile : string
        ENSEMBL ``.pep.all.fa.gz`` file in :term:`fasta` format
    outfile : string
        filename with logging information. The tablename is
        derived from ``outfile``.

    '''

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-protein_id"
        "--map=protein_id:str")

    statement = '''gunzip
    < %(infile)s
    | perl -p -e 'if ("^>") { s/ .*//};'
    | cgat fasta2fasta --method=filter
    --filter-method=min-length=1
    | cgat fasta2table --section=length
    --section=sequence
    | perl -p -e 's/id/protein_id/'
    | %(load_statement)s
    > %(outfile)s'''

    P.run()


def buildCDSFasta(infiles, outfile):
    '''output CDS sequences.

    This used to work by taking the CDNA and peptide sequence of a
    particular transcript and aligning them in order to remove any
    frameshifts.
    It relied on a deprecated library and has been removed.
    FUNCTIONALITY MISSING

    Arguments
    ---------
    infile : string
        ENSEMBL :term:`gtf` formatted file
    outfile : string
        indexed file in :term:`fasta` format with CDS sequences.

    '''
    infile_cdnas, infile_peptides_fasta = infiles

    dbname = outfile[:-len(".fasta")]

    statement = '''gunzip < %(infile_cdnas)s
    | cgat gff2fasta
        --is-gtf
        --genome=%(genome_dir)s/%(genome)s
    | cgat index_fasta
    %(dbname)s --force-output -
    > %(dbname)s.log
    '''
    P.run()


def loadGeneStats(infile, outfile):
    """compute and load gene statistics to database.

    Gene statistics are computed by :doc:`gtf2table` with the
    following counters:

    * length - gene/exon lengths
    * position - gene position
    * composition-na - gene nucleotide composition

    Parameters
    ----------
    infile : string
        A :term:`gtf` file which is output from :meth:`buildGenes`
    outfile : string
        A log file. The table name is derived from `outfile`.
        e.g. bam_stats.load
    """

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--map=gene_name:str")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2table
          --log=%(outfile)s.log
          --genome=%(genome_dir)s/%(genome)s
          --counter=position
          --counter=length
          --counter=composition-na
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


def buildExons(infile, outfile):
    '''output exons from ENSEMBL gene set.

    Remove all features from a :term:`gtf` file that are not of
    feature ``exon``.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''
    statement = '''
    gunzip < %(infile)s
    | awk '$3 == "exon"'
    | cgat gtf2gtf
    --method=remove-duplicates --duplicate-feature=gene
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


def buildCodingExons(infile, outfile):
    '''output protein coding exons from ENSEMBL gene set.

    Remove all features from a :term:`gtf` file that are not ``exon``
    and are not protein-coding.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''

    statement = '''
    zcat %(infile)s
    | cgat gtf2gtf
    --method=filter --filter-method=proteincoding
    --log=%(outfile)s.log
    | awk '$3 == "exon"'
    | cgat gtf2gtf
    --method=remove-duplicates --duplicate-feature=gene
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


def buildNonCodingExons(infile, outfile):
    '''output non-coding exons from ENSEMBL gene set.

    Remove all features from a :term:`gtf` file that are ``exon``
    and that are not protein-coding.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=filter --filter-method=proteincoding --invert-filter
    --log=%(outfile)s.log
    | awk '$3 == "exon"'
    | cgat gtf2gtf
    --method=remove-duplicates --duplicate-feature=gene
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


def buildLincRNAExons(infile, outfile):
    """output LincRNA portion of ENSEMBL geneset.

    Take all features from a :term:`gtf` file that are of feature type
    ``exon`` and that are annotated as a lincrna biotype.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    """

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=filter --filter-method=lincrna
    --log=%(outfile)s.log
    | awk '$3 == "exon"'
    | cgat gtf2gtf
    --method=remove-duplicates --duplicate-feature=gene
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


def buildCDS(infile, outfile):
    '''output CDS features from an ENSEMBL gene set.

    Take all features from a :term:`gtf` file that are of feature type
    ``CDS`` and that are annotated as protein-coding.

    Note that only the coding parts of exons are output - UTR's are
    removed.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output filename in :term:`gtf` format.

    '''
    statement = '''
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=filter --filter-method=proteincoding
    --log=%(outfile)s.log
    | awk '$3 == "CDS"'
    | cgat gtf2gtf
    --method=remove-duplicates --duplicate-feature=gene
    --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


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
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id "
        "--allow-empty-file ")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2tsv
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


def loadGeneCoordinates(infile, outfile):
    '''merge transcripts to generate the genomic coordinates per gene
    and load '''

    # TS. remove transcript_id column as this is now meaningless
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--ignore-column=transcript_id "
        "--allow-empty-file ")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2gtf
    --method=merge-transcripts
    | cgat gtf2tsv
    | %(load_statement)s
    > %(outfile)s'''

    P.run()


def loadTranscript2Gene(infile, outfile):
    '''build a map of transcript to gene from gtf file and load into database.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.
    '''
    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id ")

    statement = '''
    gunzip < %(infile)s
    | cgat gtf2tsv --output-map=transcript2gene -v 0
    | %(load_statement)s
    > %(outfile)s'''
    P.run()


def loadTranscriptStats(infile, outfile):
    '''compute and load transcript properties into database.

    The method calls :doc:`gtf2table` with the following counters:
    * length - gene/exon lengths
    * position - gene position
    * composition-na - gene nucleotide composition

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=gene_id "
        "--add-index=transcript_id "
        "--map=gene_id:str")

    statement = '''
    gunzip < %(infile)s |\
    cgat gtf2table \
          --log=%(outfile)s.log \
          --genome=%(genome_dir)s/%(genome)s \
          --reporter=transcripts \
          --counter=position \
          --counter=length \
          --counter=composition-na
    | %(load_statement)s
    > %(outfile)s'''

    P.run()


def loadProteinStats(infile, outfile):
    '''compute and load protein sequence properties into database.

    The method computes amino acid composition, length, and hash
    for each peptide sequence.

    The method calls :doc:`fasta2table` with the following counters:

    * length - protein sequence length
    * hid - protein sequence hash identifier
    * aa - protein sequence composition

    Arguments
    ---------
    infile : string
       Fiename of ENSEMBL peptide file in :term:`fasta` format.
    outfile : string
       Logfile. The table name is derived from `outfile`.

    '''

    load_statement = P.build_load_statement(
        P.toTable(outfile),
        options="--add-index=protein_id "
        "--map=protein_id:str")

    # the awk statement truncates ids ENSPXXX.1 to ENSPXXX
    # necessary for downstream compatibility (e.g. seleno list)
    statement = '''
    gunzip < %(infile)s
    | cgat fasta2fasta
    --method=filter
    --filter-method=min-length=1
    | awk 'match($0, /(>ENS[A-Z]+[0-9]+)(\.[0-9])*(.*)/, a) {print a[1], a[3]}
    !/^>/ {print}'
    | cgat fasta2table
    --log=%(outfile)s
    --sequence-type=aa
    --section=length
    --section=hid
    --section=aa
    --regex-identifier="(\S+)"
    | sed "s/^id/protein_id/"
    | %(load_statement)s
    > %(outfile)s'''

    P.run()


def buildPromotorRegions(infile, outfile, promotor_size=1000):
    '''annotate promotor regions from reference gene set.

    This method builds promotor regions for transcripts
    in an ENSEMBL gene set.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Filename in :term:`gff` format.
    promotor_size : int
       Size of the promotor region (nucleotides upstream
       of TSS).
    '''

    statement = """
    gunzip < %(infile)s
    | cgat gff2gff --method=sanitize
    --sanitize-method=genome
    --skip-missing --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | cgat gtf2gff --method=promotors
    --promotor-size=%(promotor_size)s \
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    """
    P.run()


def buildTSSRegions(infile, outfile):
    '''annotate promotor regions from reference gene set.

    This method builds promotor regions for transcripts
    in an ENSEMBL gene set.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Filename in :term:`gff` format.
    '''
    buildPromotorRegions(infile, outfile, promotor_size=1)


def buildOverlapWithEnsembl(infile, outfile, filename_bed):
    '''compute overlap of genes with intervals.

    If `filename_bed` has multiple tracks the overlap will
    be computed for each track separately.

    The output is a tab-separated table with pairs of
    overlapping features between `infile` and `filename_bed`.

    Arguments
    ---------
    infile : string
       ENSEMBL geneset in :term:`gtf` format.
    outfile : string
       Output file in :term:`tsv` format.
    filename_bed : string
       Filename in :term:`bed` format.
    '''

    statement = '''gunzip
        < %(infile)s
        | cgat gtf2gtf --method=merge-transcripts
        | cgat gff2bed --is-gtf
        | cgat bed2graph
            --output-section=name
            --log=%(outfile)s.log
            - %(filename_bed)s
        > %(outfile)s
    '''
    P.run()


def compareGeneSets(infiles, outfile):
    '''compute overlap of genes, exons and transcripts between two
    genesets.

    This method uses :mod:`scripts/diff_gtf`.

    Arguments
    ---------
    infiles : list
       Filenames of ENSEMBL genesets in :term:`gtf` format.
    outfile : string
       Output file in :term:`tsv` format.

    '''

    infiles = " ".join(infiles)
    statement = '''
        cgat diff_gtf
        %(infiles)s
    > %(outfile)s
    '''
    P.run()


def buildPseudogenes(infiles, outfile, dbhandle):
    '''build a set of pseudogenes.

    Transcripts are extracted from the GTF file and designated as
    pseudogenes if:

    * the gene_type or transcript_type contains the phrase
      "pseudo". This taken is from the database.

    * the feature is 'processed_transcript' and has similarity to
      protein coding genes. Similarity is assessed by aligning the
      transcript and peptide set against each other with exonerate_.

    Pseudogenic transcripts can overlap with protein coding
    transcripts.

    Arguments
    ---------
    infiles : list
       Filenames of ENSEMBL geneset in :term:`gtf` format
       and associated peptide sequences in :term:`fasta` format.
    outfile : filename
       Output in :term:`gtf` format with inferred or annotated
       pseudogenes.
    dbandle : object
       Database handle for extracting transcript biotypes.
    '''

    infile_gtf, infile_peptides_fasta = infiles

    # JJ - there are also 'nontranslated_CDS', but no explanation of these
    if PARAMS["genome"].startswith("dm"):
        E.warn("Ensembl dm genome annotations only contain source"
               " 'pseudogenes' - skipping exonerate step")
        statement = """zcat %(infile_gtf)s
        |awk '$2 ~ /pseudogene/'
        | gzip
        > %(outfile)s"""
        P.run()
        return

    tmpfile1 = P.getTempFilename(shared=True)

    # collect processed transcripts and save as fasta sequences
    statement = '''
    zcat %(infile_gtf)s
    | awk '$2 ~ /processed/'
    | cgat gff2fasta
            --is-gtf
            --genome-file=%(genome_dir)s/%(genome)s
            --log=%(outfile)s.log
    > %(tmpfile1)s
    '''

    P.run()

    if IOTools.isEmpty(tmpfile1):
        E.warn("no pseudogenes found")
        os.unlink(tmpfile1)
        P.touch(outfile)
        return

    model = "protein2dna"

    # map processed transcripts against peptide sequences
    statement = '''
    cat %(tmpfile1)s
    | %(cmd-farm)s --split-at-regex=\"^>(\S+)\" --chunk-size=100
    --log=%(outfile)s.log
    "exonerate --target %%STDIN%%
              --query %(infile_peptides_fasta)s
              --model %(model)s
              --bestn 1
              --score 200
              --ryo \\"%%qi\\\\t%%ti\\\\t%%s\\\\n\\"
              --showalignment no --showsugar no --showcigar no --showvulgar no
    "
    | grep -v -e "exonerate" -e "Hostname"
    | gzip > %(outfile)s.links.gz
    '''

    P.run()

    os.unlink(tmpfile1)

    inf = IOTools.openFile("%s.links.gz" % outfile)
    best_matches = {}
    for line in inf:
        peptide_id, transcript_id, score = line[:-1].split("\t")
        score = int(score)
        if transcript_id in best_matches and \
           best_matches[transcript_id][0] > score:
            continue
        best_matches[transcript_id] = (score, peptide_id)

    inf.close()

    E.info("found %i best links" % len(best_matches))
    new_pseudos = set(best_matches.keys())

    cc = dbhandle.cursor()
    known_pseudos = set([x[0] for x in cc.execute(
        """SELECT DISTINCT transcript_id
        FROM transcript_info
        WHERE transcript_biotype like '%pseudo%' OR
        gene_biotype like '%pseudo%' """)])

    E.info("pseudogenes from: processed_transcripts=%i, known_pseudos=%i, "
           "intersection=%i" % (
               (len(new_pseudos),
                len(known_pseudos),
                len(new_pseudos.intersection(known_pseudos)))))

    all_pseudos = new_pseudos.union(known_pseudos)

    c = E.Counter()

    outf = IOTools.openFile(outfile, "w")
    inf = GTF.iterator(IOTools.openFile(infile_gtf))
    for gtf in inf:
        c.input += 1
        if gtf.transcript_id not in all_pseudos:
            continue
        c.output += 1
        outf.write("%s\n" % gtf)
    outf.close()

    E.info("exons: %s" % str(c))


def buildNUMTs(infile, outfile):
    '''output set of potential nuclear mitochondrial genes (NUMTs).

    This function works by aligning the mitochondrial chromosome
    against genome using exonerate_. This can take a while.

    Arguments
    ---------
    infile : string
       Ignored.
    outfile : filename
       Output in :term:`gtf` format with potential NUMTs.

    '''
    if not PARAMS["numts_mitochrom"]:
        E.info("skipping numts creation")
        P.touch(outfile)
        return

    fasta = IndexedFasta.IndexedFasta(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"]))

    if PARAMS["numts_mitochrom"] not in fasta:
        E.warn("mitochondrial genome %s not found" % PARAMS["numts_mitochrom"])
        P.touch(outfile)
        return

    tmpfile_mito = P.getTempFilename(".")

    statement = '''
    cgat index_fasta
           --extract=%(numts_mitochrom)s
           --log=%(outfile)s.log
           %(genome_dir)s/%(genome)s
    > %(tmpfile_mito)s
    '''

    P.run()

    if IOTools.isEmpty(tmpfile_mito):
        E.warn("mitochondrial genome empty.")
        os.unlink(tmpfile_mito)
        P.touch(outfile)
        return

    format = ("qi", "qS", "qab", "qae",
              "ti", "tS", "tab", "tae",
              "s",
              "pi",
              "C")

    format = "\\\\t".join(["%%%s" % x for x in format])

    # collect all results
    min_score = 100

    statement = '''
    cat %(genome_dir)s/%(genome)s.fasta
    | %(cmd-farm)s --split-at-regex=\"^>(\S+)\" --chunk-size=1
    --log=%(outfile)s.log
    "exonerate --target %%STDIN%%
              --query %(tmpfile_mito)s
              --model affine:local
              --score %(min_score)i
              --showalignment no --showsugar no --showcigar no
              --showvulgar no
              --ryo \\"%(format)s\\n\\"
    "
    | grep -v -e "exonerate" -e "Hostname"
    | gzip > %(outfile)s.links.gz
    '''

    P.run()

    # convert to gtf
    inf = IOTools.openFile("%s.links.gz" % outfile)
    outf = IOTools.openFile(outfile, "w")

    min_score = PARAMS["numts_score"]

    c = E.Counter()

    for line in inf:
        (query_contig, query_strand, query_start, query_end,
         target_contig, target_strand, target_start, target_end,
         score, pid, alignment) = line[:-1].split("\t")

        c.input += 1
        score = int(score)
        if score < min_score:
            c.skipped += 1
            continue

        if target_strand == "-":
            target_start, target_end = target_end, target_start

        gff = GTF.Entry()
        gff.contig = target_contig
        gff.start, gff.end = int(target_start), int(target_end)
        assert gff.start < gff.end

        gff.strand = target_strand
        gff.score = int(score)
        gff.feature = "numts"
        gff.gene_id = "%s:%s-%s" % (query_contig, query_start, query_end)
        gff.transcript_id = "%s:%s-%s" % (query_contig, query_start, query_end)
        outf.write("%s\n" % str(gff))
        c.output += 1

    inf.close()
    outf.close()

    E.info("filtering numts: %s" % str(c))

    os.unlink(tmpfile_mito)


def sortGTF(infile, outfile, order="contig+gene"):
    '''sort a gtf file.

    The sorting is performed on the cluster.

    Arguments
    ---------
    infile : string
       Geneset in :term:`gtf` format.
    outfile : string
       Geneset in :term:`gtf` format.
    order : string
       Sort order. See :mod:`scripts/gtf2gtf` for valid options for
       `order`.

    '''
    if infile.endswith(".gz"):
        uncompress = "zcat"
    else:
        # wastefull
        uncompress = "cat"

    if outfile.endswith(".gz"):
        compress = "gzip"
    else:
        compress = "cat"

    job_memory = "4G"

    statement = '''%(uncompress)s %(infile)s
    | cgat gtf2gtf
    --method=sort --sort-order=%(order)s --log=%(outfile)s.log
    | %(compress)s > %(outfile)s'''

    P.run()


def buildGenomicFunctionalAnnotation(gtffile, dbh, outfiles):
    '''output a bed file with functional annotations.

    The genomic region a gene covers is taken from the `gtffile`.
    There should only be one entry per gene, i.e. exons should
    have been combined into a gene territory.

    Each entry in the output bed file is a gene territory. Bed entries
    are labeled by functional annotations associated by that gene.

    Ambiguities in territories are resolved by outputting annotations
    for all genes within a territory.

    The output file contains annotations for both GO and GOSlim. These
    are prefixed by ``go:`` and ``goslim:``.

    Arguments
    ---------
    gtffile : string
       ENSEMBL geneset in :term:`gtf` format.
    dbh : object
       Database handle to retrieve GO assignments for each gene
    outfiles : list
       Output filenames. The first is a :term:`bed` formatted file
       of gene territories. The second is a :term:`tsv` formatted
       table mapping GO terms to their description.

    '''
    outfile_bed, outfile_tsv = outfiles

    gene2region = {}
    for gtf in GTF.iterator(IOTools.openFile(gtffile, "r")):
        gid = gtf.gene_id.split(":")
        for g in gid:
            gene2region[g] = (gtf.contig, gtf.start, gtf.end, gtf.strand)

    cc = dbh.cursor()

    outf = P.getTempFile(".")
    c = E.Counter()
    term2description = {}
    for db in ('go', 'goslim'):
        for gene_id, go_id, description in cc.execute(
                "SELECT gene_id, go_id, description FROM %s_assignments" % db):
            try:
                contig, start, end, strand = gene2region[gene_id]
            except KeyError:
                c.notfound += 1
                continue
            outf.write(
                "\t".join(map(str, (
                    contig, start, end,
                    "%s:%s" % (db, go_id), 1, strand))) + "\n")
            term2description["%s:%s" % (db, go_id)] = description
    outf.close()
    tmpfname = outf.name
    statement = '''sort -k1,1 -k2,2n  < %(tmpfname)s | uniq
    | gzip > %(outfile_bed)s'''

    P.run()

    outf = IOTools.openFile(outfile_tsv, "w")
    outf.write("term\tdescription\n")
    for term, description in term2description.items():
        outf.write("%s\t%s\n" % (term, description))
    outf.close()

    os.unlink(tmpfname)


def buildGenomicContext(infiles, outfile, distance=10):
    '''build a :term:`bed` formatted file with genomic context.

    The output is a bed formatted file, annotating genomic segments
    according to whether they are any of the ENSEMBL annotations.

    The function also adds the RNA and repeats annotations from the UCSC.
    The annotations can be partially or fully overlapping.

    The annotations can be partially or fully overlapping. Adjacent
    features (less than 10 bp apart) of the same type are merged.

    Arguments
    ---------
    infiles : list
       A list of input files to generate annotations from. The contents are
       1. ``repeats``, a :term:`gff` formatted file with repeat annotations

       2. ``rna``, a :term:`gff` formatted file with small, repetetive
          RNA annotations

       3. ``annotations``, a :term:`gtf` formatted file with genomic
            annotations, see :func:`annotateGenome`.

       4. ``geneset_flat``, a flattened gene set in :term:`gtf` format, see
            :func:`buildFlatGeneSet`.

    outfile : string
       Output filename in :term:`bed` format.
    distance : int
       Merge adajcent features of the same type within this distance.

    '''

    repeats_gff, rna_gff, annotations_gtf, geneset_flat_gff, \
        cpgisland_bed, go_tsv = infiles

    tmpfile = P.getTempFilename(shared=True)
    tmpfiles = ["%s_%i" % (tmpfile, x) for x in range(6)]

    # add ENSEMBL annotations
    statement = """
    zcat %(annotations_gtf)s
    | cgat gtf2gtf
    --method=sort --sort-order=gene
    | cgat gtf2gtf
    --method=merge-exons --log=%(outfile)s.log
    | cgat gff2bed
    --set-name=gene_biotype --is-gtf
    --log=%(outfile)s.log
    | sort -k 1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_0
    """
    P.run()

    # rna
    statement = '''
    zcat %(repeats_gff)s %(rna_gff)s
    | cgat gff2bed --set-name=family --is-gtf -v 0
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_1'''
    P.run()

    # add aggregate intervals for repeats
    statement = '''
    zcat %(repeats_gff)s
    | cgat gff2bed --set-name=family --is-gtf -v 0
    | awk -v OFS="\\t" '{$4 = "repeats"; print}'
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_2'''
    P.run()

    # add aggregate intervals for rna
    statement = '''
    zcat %(rna_gff)s
    | cgat gff2bed --set-name=family --is-gtf -v 0
    | awk -v OFS="\\t" '{$4 = "repetetive_rna"; print}'
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_3 '''
    P.run()

    # add ribosomal protein coding genes
    goids = ("GO:0003735", )

    patterns = "-e %s" % ("-e ".join(goids))

    statement = '''
    zcat %(geneset_flat_gff)s
    | cgat gtf2gtf
    --map-tsv-file=<(zcat %(go_tsv)s | grep %(patterns)s | cut -f 2 | sort | uniq)
    --method=filter --filter-method=gene
    --log=%(outfile)s.log
    | cgat gff2bed
    --log=%(outfile)s.log
    | awk -v OFS="\\t" '{$4 = "ribosomal_coding"; print}'
    | sort -k1,1 -k2,2n
    | cgat bed2bed --method=merge --merge-by-name
    --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_4
    '''
    P.run()

    # CpG islands
    statement = '''
    zcat %(cpgisland_bed)s
    | awk '{printf("%%s\\t%%i\\t%%i\\tcpgisland\\n", $1,$2,$3 )}'
    > %(tmpfile)s_5
    '''
    P.run()

    # sort and merge
    # remove strand information as bedtools
    # complains if there are annotations with
    # different number of field
    files = " ".join(tmpfiles)
    statement = '''
    sort --merge -k1,1 -k2,2n %(files)s
    | cut -f 1-4
    | gzip
    > %(outfile)s
    '''
    P.run()

    for x in tmpfiles:
        os.unlink(x)
