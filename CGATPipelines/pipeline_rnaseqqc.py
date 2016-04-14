"""
====================
RNASeqQC pipeline
====================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Sailfish is utlised to rapidly estimate transcript abundance. This
requires a multi-fasta transcripts file.

For further details see http://www.cs.cmu.edu/~ckingsf/software/sailfish/

Individual tasks are enabled in the configuration file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

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

   factors=_,experiment,source,replicate,lane

Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the
apparent biases in transcript abudance within the sequence archive

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

sailfish

Code
====

ToDo
====
Documentation

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


import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineWindows as PipelineWindows
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.Pipeline as P

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
    | python %(scriptsdir)s/gtf2gtf.py
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
    | python %(scriptsdir)s/gtf2gtf.py
    --method=filter
    --filter-method=proteincoding
    --log=%(outfile)s.log
    | awk -v OFS="\\t" -v FS="\\t" '{$3="exon"; print}'
    | python %(scriptsdir)s/gtf2gtf.py
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


@follows(mkdir("nreads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"nreads.dir/\2.nreads")
def countReads(infile, outfile):
    '''Count number of reads in input files.'''
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()


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

    statement = '''python
    %(scriptsdir)s/bam2stats.py
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


@transform(buildCodingGeneSet,
           suffix(".gtf.gz"),
           ".fasta")
def buildTranscriptFasta(infile, outfile):
    """build geneset where all exons within a gene
    are merged.
    """
    dbname = outfile[:-len(".fasta")]

    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2fasta.py
    --is-gtf
    --genome=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    | python %(scriptsdir)s/index_fasta.py
    %(dbname)s --force-output -
    > %(dbname)s.log
    '''
    P.run()


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


@transform(buildCodingGeneSet,
           suffix(".gtf.gz"),
           ".tsv")
def buildTranscriptGeneMap(infile, outfile):
    """build a map of transcript ids to gene ids."""

    statement = """
    zcat %(infile)s
    |python %(scriptsdir)s/gtf2tsv.py
    --attributes-as-columns
    --output-only-attributes
    | python %(scriptsdir)s/csv_cut.py transcript_id gene_id
    > %(outfile)s"""
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
    sailfish_options += "--geneMap %s" % transcript_map

    m = PipelineMapping.Sailfish()

    statement = m.build((infile,), outfile)

    P.run()


@split(runSailfish,
       ["sailfish.dir/sailfish_transcripts.tsv.gz",
        "sailfish.dir/sailfish_genes.tsv.gz"])
def mergeSailfishResults(infiles, outfiles):

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
      if (sample_id == 1)
      {gsub(/Name/, "gene_id"); printf("sample_id\\t%%s\\n", $0); next;}}
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

    statement = '''python %(scriptsdir)s/bam2geneprofile.py
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

    infiles = [x + ".geneprofileabsolutedistancefromthreeprimeend.matrix.tsv.gz" for x in infiles]

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

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("sample_id\tfactor\tfactor_value\n")

        for sample_id, filename in enumerate(sorted(infiles)):
            sample_name, suffix = os.path.basename(filename).split(".", 1)
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
    cat %(infile)s | python %(scriptsdir)s/fasta2table.py
    --split-fasta-identifier --section=na,dn,length -v 0
    | gzip > %(outfile)s'''

    P.run()


@transform(characteriseTranscripts,
           regex("transcripts_attributes.tsv.gz"),
           add_inputs(mergeSailfishResults),
           "bias_binned_means.tsv")
def summariseBias(infiles, outfile):

    def percentage(x):
        return float(x[0])/float(x[1])

    def lin_reg_grad(x, y):
        slope, intercept, r, p, stderr = linregress(x, y)
        return slope

    attributes, genes, transcripts = infiles

    atr = pd.read_table(attributes, sep='\t', index_col="id")
    atr = atr.rename(columns={'pGC': 'GC_Content'})

    for di in itertools.product("ATCG", repeat=2):
        di = di[0]+di[1]
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
        return [(x - array_min)/(array_max-array_min) for x in array]

    def bin2floats(qcut_bin):
        qcut_bin2 = qcut_bin.replace("(", "[").replace(")", "]")
        # print "qcut_bin: ", qcut_bin, qcut_bin2
        # print type(qcut_bin)
        # print type(qcut_bin2)
        try:
            qcut_list = eval(qcut_bin2, {'__builtins__': None}, {})
            return qcut_list
        except:
            print "FAILED!!! qcut_bin: ", qcut_bin2
            return None

    def aggregate_by_factor(df, attribute, sample_names, bins, function):

        temp_dict = dict.fromkeys(sample_names, function)

        temp_dict[attribute] = function
        means_df = df[["LogTPM", "sample_id"]].groupby(
            [pd.qcut(df.ix[:, attribute], bins), "sample_id"])

        means_df = pd.DataFrame(means_df.agg(function))
        means_df.reset_index(inplace=True)

        atr_values = means_df[attribute]

        means_df.drop(attribute, axis=1, inplace=True)
        means_df["LogTPM_norm"] = norm(means_df["LogTPM"])

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


@follows(loadContextStats,
         loadBAMStats,
         loadTranscriptProfiles,
         loadSailfishResults,
         loadMetaInformation,
         loadBias)
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
