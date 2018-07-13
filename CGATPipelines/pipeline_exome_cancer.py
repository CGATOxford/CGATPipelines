
"""
======================
Exome Cancer pipeline
======================


.. todo::

   *Final filtering if SNPs/INDELs is currently done in the
   reporting. This should be handled by the pipeline. The SNP output
   would also then be passed to the mutational signature task
   *Document
   *fully make phone home/key option work - GATK public key?  Summarise
   *Indel calling (size of indels called) Example



The exome cancer pipeline imports unmapped reads from matched sample fastqs or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called and filtered


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNPs) on control samples using muTect to generate
      a "panel of normal" variants
   5a. Variant calling (SNPs) with tumour samples using muTect including
      filtering
   5b. Variant calling (indels) using Strelka
   6a. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   6b. Variant annotation with data from eBIO
   6c. Load Network of Cancer Genes (NCG) for Variant annotation in reporting


.. note::

   An optional downsampling analysis can also be performed to assess how
   coverage a control sample affects the called variants

   1. Currently the pipeline is not able to deal with replicates, i.e
      replicates will be treated seperately.



Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <patientID>-<tissue>-<replicate>.<suffix>

``patientID`` and ``tissue`` make up an :term:`experiment`, while ``replicate``
denotes the :term:`replicate` within an :term:`experiment`.
The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted
   by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

If you would like to annotate genes of interest with a particular
value in the results table, create a file call [label]_annotations.tsv
in your working directory listing all the genes. For example, to
annotate all genes identified in a previous shRNA screen, add a file
called shRNA_annoations.tsv listing the genes and the results table
will contain a column called "shRNA" with values "shRNA" and "null".

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+------------+-------------------------------------------+
|*Program*           |*Version*   |*Purpose*                                  |
+--------------------+------------+-------------------------------------------+
|Stampy              |>=0.9.0     |read mapping                               |
+--------------------+------------+-------------------------------------------+
|BWA                 |            |read mapping                               |
+--------------------+------------+-------------------------------------------+
|SAMtools            |            |filtering, SNV / indel calling             |
+--------------------+------------+-------------------------------------------+
|BEDTools            |            |filtering                                  |
+--------------------+------------+-------------------------------------------+
|sra-tools           |            |extracting reads from .sra files           |
+--------------------+------------+-------------------------------------------+
|picard              |>=1.38      |bam/sam files. The .jar files need to be in|
|                    |            |your CLASSPATH environment variable.       |
+--------------------+------------+-------------------------------------------+
|vcf-tools           |            |VCF filtering                              |
+--------------------+------------+-------------------------------------------+
|GATK                | 2.5-2      |local realignment, BQSR, variant calling   |
+--------------------+------------+-------------------------------------------+
|SNPeff              | 3.3        |                                           |
+--------------------+------------+-------------------------------------------+

Pipeline output
===============

The major output is a csvdb containing quality control information
and variant information by patientID and an html report with
similar information.

Example
=======


Code
====

"""

# load modules
from ruffus import *
# from rpy2.robjects import r as R

import numpy
import pandas as pd
import CGAT.Experiment as E
import sys
import os
import re
import shutil
import sqlite3
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineExome as PipelineExome
import CGATPipelines.PipelineBamStats as PipelineBamStats
import vcf


#########################################################################
#########################################################################


def connect():
    '''connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''
    dbh = sqlite3.connect(PARAMS["database_name"])

    return dbh


#########################################################################
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False},
    only_import=__name__ != "__main__")

PARAMS = P.PARAMS

PipelineMapping.PARAMS = PARAMS
PipelineBamStats.PARAMS = PARAMS
PipelineExome.PARAMS = PARAMS
#########################################################################

#######################################################################
# Check for design file to Match Control and Tumor input BAMs ########
#######################################################################

# This section checks for the design table and generates:
# 1. A dictionary, inputD, linking each tumour input file and its matched control,
# as specified in the design table
# 2. A pandas dataframe, df, containing the information from the
#    design table.
# 3. BAM_tumour: a list of tumour (input) bam file names following the naming guidelines
# 4. BAM_control: a list of patient matched control bam files.

# if design table is missing the input bams the df will be empty. This gets
# round the import tests

if os.path.exists("design.tsv"):
    # read the design table
    df = pd.read_csv("design.tsv", sep="\t")

    TUMOUR = list(df['BAM_tumour'].values)
    CONTROL = list(df['BAM_control'].values)

    SAMPLE_DICT = {}

    for key, value in zip(TUMOUR, CONTROL):
        SAMPLE_DICT[key] = value
else:
    E.warn("design.tsv is not located within the folder")
    TUMOUR = []
    CONTROL = []


if os.path.exists("design2.tsv"):
    # read the design table
    df = pd.read_csv("design2.tsv", sep="\t")

    TUMOUR = list(df['BAM_tumour'].values)
    CONTROL = list(df['BAM_control'].values)

    SAMPLE_DICT2 = {}

    for key, value in zip(TUMOUR, CONTROL):
        SAMPLE_DICT2[key] = value
else:
    E.warn("design2.tsv is not located within the folder")
    TUMOUR = []
    CONTROL = []

#########################################################################
# Load manual annotations
#########################################################################


@transform("*_annotations.tsv",
           suffix(".tsv"),
           ".load")
def loadManualAnnotations(infile, outfile):

    tmp = P.getTempFilename(".")

    annotation = P.snip(infile, "_annotations.tsv")

    with IOTools.openFile(tmp, "w") as outf:
        outf.write("%s\tgene_id\n" % annotation)
        with IOTools.openFile(infile, "r") as inf:
            for line in inf:
                outf.write("%s\t%s" % (annotation, line))

    P.load(tmp, outfile, options="--add-index=gene_id")
    os.unlink(tmp)

#########################################################################
# Alignment to a reference genome
#########################################################################


@follows(mkdir("bam"))
@transform(("*.fastq.1.gz", "*.fastq.gz", "*.sra"),
           regex(r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
           r"bam/\1.bam")
def mapReads(infile, outfile):
    '''Map reads to the genome using BWA, sort and index BAM file,
    generate alignment statistics and deduplicate using Picard'''

    job_threads = PARAMS["bwa_threads"]
    job_memory = PARAMS["bwa_memory"]

    if PARAMS["bwa_algorithm"] == "aln":
        m = PipelineMapping.BWA(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False)

    elif PARAMS["bwa_algorithm"] == "mem":
        m = PipelineMapping.BWAMEM(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False)
    else:
        raise ValueError("bwa algorithm '%s' not known" % algorithm)

    statement = m.build((infile,), outfile)
    print(statement)
    P.run()

#########################################################################
# Post-alignment QC
#########################################################################


@transform(mapReads,
           regex("bam/(.*).bam$"),
           r"bam/\1.dedup.bam")
def dedup(infile, outfile):
    '''Get duplicate stats from picard MarkDuplicates '''
    PipelineBamStats.buildPicardDuplicateStats(infile, outfile)


@follows(dedup)
@merge("bam/*.bam", "picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineBamStats.loadPicardDuplicateStats(infiles, outfile)

#########################################################################


@transform(dedup,
           regex("bam/(.*).bam$"),
           add_inputs(os.path.join(PARAMS["bwa_index_dir"],
                                   PARAMS["genome"])),
           r"bam/\1.picard_stats")
def buildPicardAlignStats(infiles, outfile):
    ''' build Picard alignment stats '''
    infile, reffile = infiles

    PipelineBamStats.buildPicardAlignmentStats(infile,
                                               outfile,
                                               reffile)


@follows(buildPicardAlignStats)
@merge("bam/*.picard_stats", "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineBamStats.loadPicardAlignmentStats(infiles, outfile)

#########################################################################


@transform(dedup, regex(r"bam/(\S+).dedup.bam"), r"bam/\1.dedup.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a
       bed file using Picard'''

    # TS check whether this is always required or specific to current baits
    # file

    # baits file requires modification to make picard accept it
    # this is performed before CollectHsMetrics
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]

    PipelineBamStats.buildPicardCoverageStats(
        infile, outfile, baits, regions)

    if PARAMS["zap_files"]:
        IOTools.zapFile(infile)


@follows(buildCoverageStats)
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    PipelineBamStats.loadPicardCoverageStats(infiles, outfile)

#########################################################################
#########################################################################
#########################################################################
# GATK realign bams
#########################################################################


@follows(loadCoverageStats, mkdir("gatk"))
@transform(dedup,
           regex(r"bam/(\S+).dedup.bam"),
           r"gatk/\1.readgroups.bam")
def GATKReadGroups(infile, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = re.sub(r'-\w+\.dedup.bam', '', os.path.basename(infile))
    tmpdir_gatk = P.getTempDir('.')
    job_memory = PARAMS["gatk_memory"]
    job_threads = PARAMS["gatk_threads"]
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    genome = PARAMS["bwa_indexed_genome"]

    PipelineExome.GATKReadGroups(infile, outfile, genome,
                                 library, platform,
                                 platform_unit, track)
    if PARAMS["zap_files"]:
        IOTools.zapFile(infile)


###############################################################################


@transform(GATKReadGroups,
           regex(r"gatk/(\S+).readgroups.bam"),
           r"gatk/\1.bqsr.bam")
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    #intrack = P.snip(os.path.basename(infile), ".bam")
    #outtrack = P.snip(os.path.basename(outfile), ".bam")
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["bwa_indexed_genome"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]

    PipelineExome.GATKBaseRecal(infile, outfile, genome, intervals,
                                padding, dbsnp, solid_options)
    if PARAMS["zap_files"]:
        IOTools.zapFile(infile)


#########################################################################
#########################################################################
#########################################################################
# Variant Calling
#########################################################################


@follows(mkdir("normal_panel_variants", GATKBaseRecal))
@transform(GATKBaseRecal,
           regex(r"gatk/(\S+)-%s-(\S+).bqsr.bam" %
                 PARAMS["sample_control"]),
           r"normal_panel_variants/\1_normal_mutect.vcf")
def callControlVariants(infile, outfile):
    '''run mutect to call snps in control sample'''

    basename = P.snip(outfile, "_normal_mutect.vcf")
    mutect_log = basename + ".log"

    cosmic = PARAMS["mutect_cosmic"]
    dbsnp = PARAMS["mutect_dbsnp"]
    roi_intervals = PARAMS["roi_intervals"]
    job_threads = PARAMS['mutect_threads']
    job_memory = PARAMS['mutect_memory']

    genome = PARAMS["bwa_indexed_genome"]

    PipelineExome.MuTect2Caller(outfile, mutect_log, genome,
                                cosmic, dbsnp, job_memory, job_threads,
                                roi_intervals, infile)


@transform(callControlVariants,
           regex(r"normal_panel_variants/(\S+).vcf"),
           r"normal_panel_variants/\1_slim.vcf.gz")
def indexControlVariants(infile, outfile):
    '''index control vcf for intersection by vcftools'''
    '''tabix is a tool that allows to perform fast interval queries based on tab delimited interval file'''

    outfile = P.snip(outfile, ".gz")

    statement = '''cut -f1-8 %(infile)s > %(outfile)s;
                   bgzip -f %(outfile)s;
                   tabix -f %(outfile)s.gz'''
    P.run()


# paramaterise vcf intersection (number of req. observations - currently 1)
@merge(indexControlVariants,
       "normal_panel_variants/combined.vcf")
def mergeControlVariants(infiles, outfile):
    ''' intersect control vcfs to generate a panel of normals for mutect'''
    '''outputs positions present in at least one file'''
    infiles = " ".join(infiles)

    statement = '''module load vcftools/0.1.14;
                   module load perl;
                   export PERL5LIB=/package/vcftools/0.1.14/lib/site_perl/5.18.1;
                   vcf-isec -n +1 %(infiles)s
                   > %(outfile)s'''
    P.run()


@follows(GATKBaseRecal)
@collate(GATKBaseRecal,
         regex(r"gatk/(CM[0-9]{4})(\S+).bqsr.bam"),
         r"gatk/\1.pt")
def patientID(infiles, outfile):
    '''makes and empty file for patient ID'''
    '''patient sample names should start with CM'''
    '''might need to change it for different patient names'''

    to_cluster = False
    statement = '''touch %(outfile)s'''

    P.run()


@follows(mkdir("variants"), patientID)
@transform(patientID,
           regex(r"gatk/(.*).pt"),
           r"variants/\1.mutect.vcf")
def runMutect(infile, outfile):
    '''calls somatic SNPs and indels using MuTect2'''

    #base = P.snip(os.path.basename(infile), ".pt")
    infile_tumour = P.snip(infile, ".pt") + "-Tumour-R1.bqsr.bam"
    #infile_tumour_key = P.snip(os.path.basename(infile), ".pt") + "-Tumour-R1"
    #infile_control = "gatk/" + SAMPLE_DICT[infile_tumour_key] + ".bqsr.bam"
    infile_control = P.snip(
        infile_tumour, "-Tumour-R1.bqsr.bam") + "-Control-R1.bqsr.bam"

    basename = P.snip(outfile, ".mutect.vcf")
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, minReadsPerAlignmentStart, tumor_LOD, normal_LOD, strand_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_minreadsperalignmentstart"],
         PARAMS["mutect_lod"], PARAMS["mutect_normal_lod"], PARAMS["mutect_strand_lod"])

    job_threads = PARAMS['mutect_threads']
    job_memory = PARAMS['mutect_memory']

    roi_intervals = PARAMS["roi_intervals"]
    genome = PARAMS["bwa_indexed_genome"]

    PipelineExome.MuTect2Caller(
        outfile, mutect_log, genome,
        cosmic, dbsnp,
        job_memory, job_threads, roi_intervals,
        infile_tumour=infile_tumour, infile_control=infile_control,
        quality=quality, max_alt_qual=max_alt_qual,
        max_alt=max_alt, max_fraction=max_fraction,
        minReadsPerAlignmentStart=minReadsPerAlignmentStart,
        tumor_LOD=tumor_LOD, normal_LOD=normal_LOD, strand_LOD=strand_LOD)


##########################################################################
##########################################################################
##########################################################################
# repeat mutect in reverse and on subsampled control bam as quality control
##########################################################################
# this analysis should be part of an optional check of mutect parameters
# mutect paramters should be identical to the runMutect function above
# splitMergedRealigned replaced by GATKBaseRecal


@follows(patientID)
@transform(patientID,
           regex(r"gatk/(.*).pt"),
           r"variants/\1.mutect.reverse.vcf")
def runMutectReverse(infile, outfile):
    '''calls somatic SNPs and indels using MuTect2'''

    #base = P.snip(os.path.basename(infile), ".pt")
    infile_control = P.snip(infile, ".pt") + "-Tumour-R1.bqsr.bam"
    infile_tumour_key = P.snip(os.path.basename(infile), ".pt") + "-Tumour-R1"
    infile_tumour = "gatk/" + SAMPLE_DICT[infile_tumour_key] + ".bqsr.bam"

    basename = P.snip(outfile, ".mutect.reverse.vcf")
    mutect_log = basename + ".reverse.log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, strand_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_lod"], PARAMS["mutect_strand_lod"])

    job_threads = PARAMS['mutect_threads']
    job_memory = PARAMS['mutect_memory']

    roi_intervals = PARAMS["roi_intervals"]
    genome = PARAMS["bwa_indexed_genome"]

    PipelineExome.MuTect2Caller(
        outfile,
        mutect_log,
        genome,
        cosmic,
        dbsnp,
        job_memory,
        job_threads,
        roi_intervals,
        infile_tumour=infile_tumour,
        infile_control=infile_control,
        quality=quality,
        max_alt_qual=max_alt_qual,
        max_alt=max_alt,
        max_fraction=max_fraction,
        tumor_LOD=tumor_LOD,
        strand_LOD=strand_LOD)


'''@follows(mergeControlVariants)
@transform(GATKBaseRecal,
           regex(r"bam/(\S+)-%s-(\S+).realigned.split.bqsr.bam" %
                 PARAMS["sample_control"]),
           add_inputs(mergeControlVariants),
           r"variants/\1.mutect.reverse.snp.vcf")
def runMutectReverse(infiles, outfile):'''
'''Use control as tumor and vis versa to estimate false positive rate'''
'''    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        PARAMS["sample_control"], PARAMS["sample_tumour"])

    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    basename = P.snip(outfile, ".mutect.reverse.snp.vcf")
    call_stats_out = basename + "_call_stats.reverse.out"
    coverage_wig_out = basename + "_coverage.reverse.wig"
    mutect_log = basename + ".reverse.log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineExome.mutectSNPCaller(infile, outfile, mutect_log, genome,
                                  cosmic, dbsnp, call_stats_out,
                                  PARAMS['mutect_memory'],
                                  PARAMS['mutect_threads'],
                                  quality, max_alt_qual,
                                  max_alt, max_fraction, tumor_LOD,
                                  normal_panel, infile_tumour)'''


# generalise the functions below
# 1. identify sample with highest coverage in control
# - should this check coverage in tumour also?
# 2. subset control bam
# 3. run mutect calling function with subset against unsubsetted tumour
# 4. summary table

#adeno_bam = "bam/NU16C-Control-1.realigned.bqsr.bam"
#adeno_bam = "gatk/CM0106-Control-R1.bqsr.bam"
# GATKBaseRecal is used instead of adeno_bam

@subdivide(GATKBaseRecal,
           regex("(\S+).bqsr.bam"),
           [r"\1.0.1.subd.bqsr.bam",
            r"\1.0.2.subd.bqsr.bam",
            r"\1.0.3.subd.bqsr.bam",
            r"\1.0.4.subd.bqsr.bam",
            r"\1.0.5.subd.bqsr.bam",
            r"\1.0.6.subd.bqsr.bam",
            r"\1.0.7.subd.bqsr.bam",
            r"\1.0.8.subd.bqsr.bam",
            r"\1.0.9.subd.bqsr.bam",
            r"\1.1.0.subd.bqsr.bam"])
def subsetBqsrBam(infile, outfiles):
    statements = []
    n = 0
    for fraction in numpy.arange(0.1, 1.1, 0.1):
        outfile = outfiles[n]
        n += 1
        statement = '''samtools view -s %(fraction)s -b %(infile)s
                     > %(outfile)s 2> %(outfile)s.stderr'''
        P.run()


@transform(subsetBqsrBam,
           suffix(".subd.bqsr.bam"),
           ".subd.bqsr.bam.bai")
def indexSubsets(infile, outfile):
    sorted_bam = "gatk/sorted." + os.path.basename(infile)
    statement = '''samtools sort -T ${TMPDIR}/${USER} -o %(sorted_bam)s %(infile)s ;
                   samtools index %(sorted_bam)s'''
    P.run()


@collate(subsetBqsrBam,
         regex(r"gatk/(CM[0-9]{4})(\S+)-R1(\S+).subd.bqsr.bam"),
         r"gatk/\1-\3.pt")
def patientID2(infiles, outfile):
    '''makes and empty file for patient ID and subset sample'''
    '''patient sample names should start with CM'''
    '''might need to change it for different patient names'''

    to_cluster = False
    statement = '''touch %(outfile)s'''

    P.run()


@transform(patientID2,
           regex(r"gatk/(\S+)-.(\S+).pt"),
           r"variants/\1-.\2-downsampled.mutect.vcf")
def runMutectOnDownsampled(infile, outfile):
    '''calls somatic SNPs and indels using MuTect2'''

    #base = P.snip(os.path.basename(infile), ".pt")
    m = re.search(r'gatk/(\S+)-.(\S+).pt', infile)
    m.group(1)
    m.group(2)
    infile_tumour = "gatk/" + "sorted." + \
        m.group(1) + "-Tumour-R1." + m.group(2) + ".subd.bqsr.bam"
    infile_tumour_key = "sorted." + \
        m.group(1) + "-Tumour-R1." + m.group(2) + ".subd.bqsr.bam"
    infile_control = "gatk/" + SAMPLE_DICT2[infile_tumour_key]

    basename = P.snip(outfile, "-downsampled.mutect.vcf")
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, strand_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_lod"], PARAMS["mutect_strand_lod"])

    job_threads = PARAMS['mutect_threads']
    job_memory = PARAMS['mutect_memory']

    roi_intervals = PARAMS["roi_intervals"]
    genome = PARAMS["bwa_indexed_genome"]

    PipelineExome.MuTect2Caller(
        outfile,
        mutect_log,
        genome,
        cosmic,
        dbsnp,
        job_memory,
        job_threads,
        roi_intervals,
        infile_tumour=infile_tumour,
        infile_control=infile_control,
        quality=quality,
        max_alt_qual=max_alt_qual,
        max_alt=max_alt,
        max_fraction=max_fraction,
        tumor_LOD=tumor_LOD,
        strand_LOD=strand_LOD)


'''@follows(indexSubsets)
@transform(subsetBqsrBam,
           regex(r"bam/(\S+)-%s-1.realigned.(\S+).bqsr.bam" %
                 PARAMS["sample_control"]),
           add_inputs(mergeControlVariants),
           r"variants/\1-downsampled-\2.mutect.snp.vcf")'''
'''def runMutectOnDownsampled(infiles, outfile):'''
'''call somatic SNPs using MuTect on downsampled bams'''
'''    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        PARAMS["sample_control"], PARAMS["sample_tumour"])
    basename = P.snip(outfile, "_normal_mutect.vcf")

    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineExome.mutectSNPCaller(infile_tumour, outfile, mutect_log, genome,
                                  cosmic, dbsnp, call_stats_out,
                                  PARAMS['mutect_memory'], PARAMS[
                                      'mutect_threads'],
                                  quality, max_alt_qual,
                                  max_alt, max_fraction, tumor_LOD,
                                  normal_panel, infile)'''

##############################################################################
##############################################################################
##############################################################################
# Variant Annotation
##############################################################################


'''@collate(GATKBaseRecal,
         regex(r"gatk/(\S+)-%s-(\S+).bqsr.bam" %
                 PARAMS["sample_control"]),
         r"gatk/\1.list")
def listOfBAMs(infiles, outfile):'''
'''generates a file containing a list of BAMs for each patient,
       for use in variant calling'''
'''    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            infile_tumour = infile.replace(
                PARAMS["sample_control"], PARAMS["sample_tumour"])
            outf.write(infile + '\n')
            outf.write(infile_tumour + '\n')'''


@follows(runMutect, mkdir("annotations"))
@transform(runMutect,
           regex(r"variants/(\S+).mutect.vcf"),
           r"annotations/\1.mutect.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    job_memory = "4G"
    job_threads = 2

    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''snpEff -Xmx4G eff
                   -c %(config)s -v %(snpeff_genome)s -o gatk
                   %(infile)s > %(outfile)s'''
    P.run()


#########################################################################
# Annotate SNP and INDEL variants GATK
#########################################################################

# Need to check whether variant annotator is using both bams
# from a single patient?
# should just be the tumour bam or else scores will be wrong!

@follows(annotateVariantsSNPeff)
@transform(patientID,
           regex(r"gatk/(.*).pt"),
           r"annotations/\1.mutect.annotated.vcf")
def variantAnnotator(infile, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''

    vcffile = "variants/" + \
        P.snip(os.path.basename(infile), ".pt") + ".mutect.vcf"
    bamfile = "gatk/" + P.snip(os.path.basename(infile),
                               ".pt") + "-Tumour-R1.bqsr.bam"
    snpeff_file = "annotations/" + \
        P.snip(os.path.basename(infile), ".pt") + ".mutect.snpeff.vcf"

    dbsnp = PARAMS["gatk_dbsnp"]
    genome = PARAMS["bwa_indexed_genome"]
    annotations = PARAMS["annotation_variant_annotations"]
    PipelineExome.variantAnnotator(vcffile, bamfile, outfile, genome,
                                   dbsnp, annotations, snpeff_file)


@transform(variantAnnotator,
           regex(r"annotations/(\S+).mutect.annotated.vcf"),
           r"annotations/\1.mutect.annotated.cosmic.vcf")
def annotateVariantsCosmic(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = PARAMS["annotation_memory"]
    job_threads = PARAMS["annotation_threads"]

    cosmic = PARAMS["annotation_cosmic"]
    statement = '''SnpSift annotate
                %(cosmic)s
                %(infile)s > %(outfile)s;'''

    P.run()

##################################################################


@transform(annotateVariantsCosmic,
           regex(r"annotations/(\S+).mutect.annotated.cosmic.vcf"),
           r"annotations/\1.dbnsfp.snpsift.vcf")
def annotateVariantsDBNSFP(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = PARAMS["annotation_memory"]
    job_threads = PARAMS["annotation_threads"]
    dbNSFP = PARAMS["annotation_dbnsfp"]
    config = PARAMS["annotation_config"]
    genome = PARAMS["annotation_snpsift_genome"]
    dbN_annotators = PARAMS["annotation_dbnsfpannotators"]
    if len(dbN_annotators) == 0:
        annostring = ""
    else:
        annostring = "-f %s" % dbN_annotators

    statement = """SnpSift dbnsfp -db %(dbNSFP)s -config %(config)s  -v %(infile)s
                   %(annostring)s >
                   %(outfile)s;"""
    P.run()

##################################################################


@transform(annotateVariantsDBNSFP,
           regex(r"annotations/(\S+).dbnsfp.snpsift.vcf"),
           r"annotations/\1.exac.snpsift.vcf")
def annotateVariantsExAC(infile, outfile):
    '''Add annotations using SNPsift, Exome Aggregation Consortium (ExAC)'''
    job_memory = PARAMS["annotation_memory"]
    job_threads = PARAMS["annotation_threads"]
    exac = PARAMS["annotation_exac"]
    statement = """SnpSift annotate
                %(exac)s
                %(infile)s > %(outfile)s;"""
    P.run()

####################################################################


@transform(annotateVariantsExAC,
           regex(r"annotations/(\S+).exac.snpsift.vcf"),
           r"annotations/\1.1000G.snpsift.vcf")
def annotateVariants1000G(infile, outfile):
    '''Add annotations using SNPsift'''
    job_memory = PARAMS["annotation_memory"]
    job_threads = PARAMS["annotation_threads"]

    vcfs = []
    for f in os.listdir(PARAMS["annotation_tgdir"]):
        if f.endswith(".vcf.gz"):
            if PARAMS['test'] == 1:
                if "chr14" in f:
                    vcfs.append("%s/%s" % (PARAMS['annotation_tgdir'], f))
            else:
                vcfs.append("%s/%s" % (PARAMS['annotation_tgdir'], f))

    T = P.getTempFilename(".")
    shutil.copy(infile, T)
    tempin = T
    tempout = P.getTempFilename(".")

    for vcf in vcfs:
        statement = """SnpSift annotate
                       %(vcf)s
                       %(tempin)s 1> %(tempout)s 2>> %(tempout)s.stderr;
                       mv %(tempout)s %(tempin)s"""
        P.run()

    shutil.move(tempin, outfile)

################################################################


@transform(annotateVariants1000G,
           regex(r"annotations/(\S+).1000G.snpsift.vcf"),
           r"annotations/\1.vep.vcf")
def annotateVariantsVEP(infile, outfile):
    '''
    Adds annotations as specified in the pipeline.ini using Ensembl
    variant effect predictor (VEP).
    '''
    # infile - VCF
    # outfile - VCF with vep annotations
    job_memory = PARAMS["annotation_memory"]
    job_threads = 4

    VEP = PARAMS["annotation_vepannotators"].split(",")
    vep_annotators = PARAMS["annotation_vepannotators"]
    vep_path = PARAMS["annotation_veppath"]
    vep_cache = PARAMS["annotation_vepcache"]
    vep_species = PARAMS["annotation_vepspecies"]
    vep_assembly = PARAMS["annotation_vepassembly"]
    if len(vep_annotators) != 0:
        annostring = vep_annotators
        statement = '''vep --cache --dir %(vep_cache)s --vcf
                       --species %(vep_species)s
                       --fork 2
                       --assembly %(vep_assembly)s --input_file %(infile)s
                       --output_file %(outfile)s --force_overwrite
                       %(annostring)s --offline;'''
        P.run()

##############################################################################
##############################################################################
##############################################################################
############################################################################
# split the table annotations into columns


@follows(mkdir("variant_tables"))
@transform(annotateVariantsVEP,
           regex(r"annotations/(\S+).vep.vcf"),
           r"variant_tables/\1_TUMOR.tsv")
def makeAnnotationsTablesTumor(infile, outfile):
    '''
    Converts the vcf generated with MuTect2 into a table containing only
    positions called as variants in that sample and INFO field split into columns.
    '''

    #TF = P.getTempFilename(".")
    samplename = "TUMOR"

    vcf_reader = vcf.Reader(open(infile))
    # requires installation of pyvcf

    list_IDs = []
    list_desc = []
    list_cstring = []
    for k, v in vcf_reader.formats.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    for k, v in vcf_reader.infos.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    out = open(outfile, "w")
    out.write('''POS\tCHROM\tQUAL\tID\tFILTER\tREF1\tALT\tGT\t%s\nposition\
              \tchromosome\tquality\tid\tfilter\tref\
              \talt\tgenotype\t%s\n''' % ("\t".join(list_IDs), "\t".join(list_desc)))
    out.close()
    cstring = "\\t".join(list_cstring)
    cstring = "%POS\\t%CHROM\\t%QUAL\\t%ID\\t%FILTER\\t%REF\\t\
               %ALT\\t[%GT]\\t" + cstring

    if PARAMS['test'] == 1:
        statement = '''bcftools query
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    else:
        statement = '''bcftools query -s %(samplename)s
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    P.run()


############################################################################


@transform(annotateVariantsVEP,
           regex(r"annotations/(\S+).vep.vcf"),
           r"variant_tables/\1_NORMAL.tsv")
def makeAnnotationsTablesNormal(infile, outfile):
    '''
    Converts the vcf generated with MuTect2 into a table containing only
    positions called as variants in that sample and INFO field split into columns.
    '''

    #TF = P.getTempFilename(".")
    samplename = "NORMAL"

    vcf_reader = vcf.Reader(open(infile))
    # requires installation of pyvcf

    list_IDs = []
    list_desc = []
    list_cstring = []
    for k, v in vcf_reader.formats.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    for k, v in vcf_reader.infos.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    out = open(outfile, "w")
    out.write('''POS\tGT\tAB\tAD\tAF\tALT_F1R2\tALT_F2R1\tDP\tFOXOG\tPGT\tPID\tQSS\tREF_F1R2\
              \tREF_F2R1\nposition\tNormal_genotype\
                 \tNormal_Allele_balance_for_each_het_genotype\
                 \tNormal_Allelic_depths_for_the_ref_and_alt_alleles_in_the_order_listed\
                 \tNormal_Allele_fraction_of_the_event_in_the_tumor\
                 \tNormal_Count_of_reads_in_F1R2_pair_orientation_supporting_the_alternate_allele\
                 \tNormal_Count_of_reads_in_F2R1_pair_orientation_supporting_the_alternate_allele\
                 \tNormal_Approximate_read_depth_(reads_with_MQ=255_or_with_bad_mates_are_filtered)\
                 \tNormal_Fraction_of_alt_reads_indicating_OxoG_error\
                 \tNormal_Physical_phasing_haplotype_information,_describing_how_the_alternate_alleles_are_phased_in_relation_to_one_another\
                 \tNormal_Physical_phasing_ID_information,_where_each_unique_ID_within_a_given_sample_(but_not_across_samples)_connects_records_within_a_phasing_group\
                 \tNormal_Sum_of_base_quality_scores_for_each_allele\
                 \tNormal_Count_of_reads_in_F1R2_pair_orientation_supporting_the_reference_allele\
                 \tNormal_Count_of_reads_in_F2R1_pair_orientation_supporting_the_reference_allele\n''')
    out.close()
    cstring = "\\t".join(list_cstring)
    cstring = "%POS\\t[%GT]\\t[%AB]\\t[%AD]\\t[%AF]\\t[%ALT_F1R2]\\t[%ALT_F2R1]\\t[%DP]\\t[%FOXOG]\\t[%PGT]\\t[%PID]\\t[%QSS]\\t[%REF_F1R2]\\t[%REF_F2R1]"

    if PARAMS['test'] == 1:
        statement = '''bcftools query
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    else:
        statement = '''bcftools query -s %(samplename)s
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    P.run()


########################################################################


@transform(annotateVariantsVEP,
           regex(r"annotations/(\S+).vep.vcf"),
           r"variant_tables/\1_NORMAL_Sel.tsv")
def makeAnnotationsTablesNormalSel(infile, outfile):
    '''
    Converts the vcf generated with MuTect2 into a table containing only
    positions called as variants in that sample and INFO field split into columns.
    '''

    #TF = P.getTempFilename(".")
    samplename = "NORMAL"

    vcf_reader = vcf.Reader(open(infile))
    # requires installation of pyvcf

    list_IDs = []
    list_desc = []
    list_cstring = []
    for k, v in vcf_reader.formats.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    for k, v in vcf_reader.infos.items():
        if k != "Samples":
            list_IDs.append(k)
            list_cstring.append("[%%%s]" % k)
            list_desc.append(v.desc.strip().replace(" ", "_"))

    out = open(outfile, "w")
    out.write('''POS\tCHROM\tGT\tAB\tAD\nposition\tchromosome\tNormal_genotype\
              \tNormal_Allele_balance_for_each_het_genotype\
              \tNormal_Allelic_depths_for_the_ref_and_alt_alleles_in_the_order_listed\n''')

    out.close()
    cstring = "\\t".join(list_cstring)
    cstring = "%POS\\t%CHROM\\t[%GT]\\t[%AB]\\t[%AD]"

    if PARAMS['test'] == 1:
        statement = '''bcftools query
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    else:
        statement = '''bcftools query -s %(samplename)s
                   -f '%(cstring)s\\n'
                   %(infile)s >> %(outfile)s'''
    P.run()

##########################################################################


@follows(makeAnnotationsTablesNormal)
@transform(makeAnnotationsTablesTumor,
           regex("variant_tables/(\S+)_TUMOR.tsv"),
           r"variant_tables/\1.combined.tsv")
def mergeTables(infile, outfile):
    ''' merge TUMOR and CONTROL tables '''

    job_memory = PARAMS["annotation_memory"]
    job_threads = 4

    infile_TUMOR = infile
    infile_NORMAL = P.snip(infile_TUMOR, "_TUMOR.tsv") + "_NORMAL_Sel.tsv"

    statement = '''join -j 1 -t $'\t' %(infile_NORMAL)s %(infile_TUMOR)s > %(outfile)s'''

    P.run()

#############################################################################


@transform(variantAnnotator,
           regex("annotations/(\S+).mutect.annotated.vcf"),
           r"annotations/\1.mutect.annotated.filtered.vcf")
def filterMutect(infile, outfile):
    ''' filter MuTect2 snps and indels using allele frequencies '''

    logfile = outfile.replace(".vcf", ".log")

    min_t_alt = PARAMS["filter_minimum_tumor_allele"]
    min_t_alt_freq = PARAMS["filter_minimum_tumor_allele_frequency"]
    min_n_depth = PARAMS["filter_minimum_normal_depth"]
    max_n_alt_freq = PARAMS["filter_maximum_normal_allele_frequency"]
    min_ratio = PARAMS["filter_minimum_ratio"]

    PipelineExome.filterMutect(
        infile, outfile, logfile,
        min_t_alt, min_n_depth, max_n_alt_freq,
        min_t_alt_freq, min_ratio)

##############################################################################
# Intersect filtered SNPs and INDELs
##############################################################################

# also included filterIndels


@mkdir("intersection.dir")
@collate((filterMutect),
         regex(r"variants/(\S+)\.(\S+).annotated.filtered.vcf"),
         r"intersection.dir/overlap_\2_heatmap.png")
def intersectHeatmap(infiles, outfile):
    ''' intersect DE test_ids across the different quantifiers'''

    PipelineExome.intersectionHeatmap(infiles, outfile)

#########################################################################
#########################################################################
# convert vcf to tsv files and load into database


@transform(filterMutect,
           regex("variants/(\S+).annotated.filtered.vcf"),
           r"variants/\1.annotated.filtered.tsv")
def snpvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    statement = '''GenomeAnalysisTK
                   -T VariantsToTable -R %(bwa_index_dir)s/%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -GF GT -GF AD -GF SS -GF FA -GF AB -GF DP
                   -o %(outfile)s'''
    P.run()


# was filterIndels
@transform(filterMutect,
           regex("variants/(\S+).annotated.filtered.vcf"),
           r"variants/\1.annotated.filtered.tsv")
def indelvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    statement = '''GenomeAnalysisTK
                   -T VariantsToTable -R %(bwa_index_dir)s/%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -F TQSI -F TSQI_NT -F DP -F IC -F IHP -F NT
                   -F QSI -F QSI_NT -F RC -F RU -F SGT
                   -GF DP -GF DP2 -GF DP50 -GF SUBDP50 -GF TAR -GF TIR -GF TOR
                   -o %(outfile)s'''
# The parameter allowMissingData is deprecated.
    P.run()


@transform([snpvcfToTable,
            indelvcfToTable],
           regex(r"variants/(\S+).annotated.filtered.tsv"),
           r"variants/\1_annotated.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''

    if infile.endswith("indels.annotated.filtered.tsv"):
        indices = "CHROM,POS,SNPEFF_GENE_NAME"
    elif infile.endswith("mutect.snp.annotated.filtered.tsv"):
        indices = "CHROM,POS,SNPEFF_GENE_NAME"

    P.load(infile, outfile, options="--add-index=%(indices)s" % locals())


#########################################################################
#########################################################################
#########################################################################
# vcf statistics -   this only summarises the nucleotide changes
# this currently does not provide useful output!

# task included variantAnnotatorIndels


@transform((variantAnnotator),
           regex(r"variants/(\S+).vcf"),
           r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    statement = '''vcf-stats %(infile)s
                   > %(outfile)s 2>>%(outfile)s.log;'''
    P.run()


@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    csv2db_options = PARAMS["csv2db_options"]
    E.info("Loading vcf stats...")
    statement = '''cgat vcfstats2db
                   %(filenames)s >> %(outfile)s; '''
    statement += '''cat vcfstats.txt |
                    cgat csv2db %(csv2db_options)s
                    --allow-empty-file --add-index=track --table=vcf_stats
                    >> %(outfile)s; '''
    P.run()

#########################################################################


@transform(runMutect,
           suffix(".mutect.snp.vcf"),
           "_mutect_filtering_summary.tsv")
def summariseFiltering(infile, outfile):
    infile = infile.replace(".mutect.snp.vcf", "_call_stats.out")

    PipelineExome.parseMutectCallStats(infile, outfile, submit=True)


@transform(summariseFiltering,
           regex(r"variants/(\S+)_mutect_filtering_summary.tsv"),
           r"variants/\1_mutect_filtering_summary.load")
def loadMutectFilteringSummary(infile, outfile):
    '''Load mutect extended output into database'''

    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   cgat csv2db
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@originate("eBio_studies.tsv")
def defineEBioStudies(outfile):
    ''' For the cancer types specified in pipeline.ini, identify the
    relevent studies in eBio '''

    cancer_types = PARAMS["annotation_ebio_cancer_types"]

    PipelineExome.defineEBioStudies(cancer_types, outfile, submit=False)

# also included filterIndels


@transform(defineEBioStudies,
           suffix("eBio_studies.tsv"),
           add_inputs(filterMutect),
           "eBio_studies_gene_frequencies.tsv")
def extractEBioinfo(infiles, outfile):
    '''find the number of mutations identitified in previous studies (ebio_ids)
    for the mutated genes in the annotated vcfs'''

    eBio_ids = infiles[0]
    vcfs = infiles[1:]

    PipelineExome.extractEBioinfo(eBio_ids, vcfs, outfile, submit=False)


@transform(extractEBioinfo,
           suffix(".tsv"),
           ".load")
def loadEBioInfo(infile, outfile):
    '''load the frequencies from the eBIO portal'''

    P.load(infile, outfile, options="--add-index=gene")

#########################################################################
#########################################################################
#########################################################################
# load Network of Cancer Genes table

# parameterise file location:


@originate("cancergenes.load")
def loadNCG(outfile):
    '''Load NCG into database'''

    infile = PARAMS["cancergenes_table"]
    # infile = "/ifs/projects/proj053/backup/NCG/cancergenes2016.tsv"

    P.load(infile, outfile, options="--add-index=symbol")

#########################################################################
#########################################################################
#########################################################################
# analyse mutational siganture of filtered variants


@merge(filterMutect,
       ["variants/mutational_signature.tsv",
        "variants/mutational_signature_table.tsv"])
def mutationalSignature(infiles, outfiles):

    PipelineExome.compileMutationalSignature(
        infiles, outfiles)


@transform(mutationalSignature,
           suffix(".tsv"),
           ".load")
def loadMutationalSignature(infiles, outfile):
    outfile2 = re.sub(".load", "_table.load", outfile)
    P.load(infiles[0], outfile)
    P.load(infiles[1], outfile2)


#########################################################################
#########################################################################
#########################################################################


@follows(defineEBioStudies)
def test():
    pass


@follows(runMutectOnDownsampled,
         runMutectReverse)
def TestMutect():
    '''This target runs function which can be used to assess the chosen
    mutect parameters'''


# @follows(loadROI,
#         loadROI2Gene)
# def loadMetadata():
#    pass


@follows(mapReads)
def mapping():
    pass


@follows(dedup,
         loadPicardDuplicateStats,
         buildPicardAlignStats,
         loadPicardAlignStats,
         buildCoverageStats,
         loadCoverageStats)
def postMappingQC():
    pass


@follows(GATKReadGroups,
         GATKBaseRecal)
def gatk():
    pass


@follows(patientID,
         runMutect)
def callVariants():
    pass


@follows(loadVariantAnnotation)
def tabulation():
    pass


@follows(annotateVariantsSNPeff,
         variantAnnotator,
         annotateVariantsCosmic,
         annotateVariantsDBNSFP,
         annotateVariantsExAC,
         annotateVariants1000G,
         annotateVariantsVEP)
def annotation():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(postMappingQC,
         gatk,
         callVariants,
         loadManualAnnotations,
         annotation,
         loadMutectFilteringSummary,
         loadVariantAnnotation,
         loadCoverageStats,
         loadPicardAlignStats,
         loadNCG,
         loadMutationalSignature,
         loadEBioInfo,
         intersectHeatmap)
def full():
    pass

#########################################################################
#########################################################################
#########################################################################


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


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
