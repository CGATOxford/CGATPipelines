"""
====================
Exome pipeline
====================

:Author: David Sims & Kath Fawcett
:Release: $Id$
:Date: |today|
:Tags: Python

The exome pipeline imports unmapped reads from one or more fastq or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called, annotated,
and filtered according to various inheritance models (de novo,
dominant, and recessive).


   1. Align to genome using gapped alignment (BWA-MEM)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK and deduplication in Picard
   4. Calculate the ratio of reads on the X and Y chromosomes to assert sex
   5. Variant calling in families (SNVs & indels) using GATK HaplotypeCaller
   6. HapMap genotyping in individuals using GATK HaplotypeCaller
   7. Comparison of Hapmap genotypes to assess relatedness (VCFtools)
   8. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   9. Variant quality score recalibration (GATK)
   10. Flags variants within genes of interest (such as known disease genes)
      (GATK) (optional)
   11. Filters potential de novo variants
   12. Filters potential de novo variants using lower stringency criteria
   13. Filters potential dominant mutations
   14. Filters potential homozygous recessive mutations
   15. Filters potential compound heterozygous mutations
   16. Generates summary statistics for unfiltered vcf file
   17. Generates report

.. note::

   1. Great care should be taken in interpreting lower stringency de
      novo variants as it is expected almost all will be false
      positives. Users should examine them manually in the csvdb
      database and html report.

   2. Great care should be taken when interpreting compound
      heterozygous changes.  Gemini is very permissive and users
      should examine the genotypes in the csvdb table to make sure
      they are consistent with a recessive inheritance pattern.

To do:
   1. Allow users to add other training sets for variant quality score
      recalibration
   2. Allow users to add annotations using SnpSift

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

   <family>_<sample>-<condition>-<replicate>.<suffix>

``family`` = "Single", "Trio", or "Multiplex" followed by numerical
identifier.  ``sample`` and ``condition`` make up an
:term:`experiment`, while ``replicate`` denotes the :term:`replicate`
within an :term:`experiment`.  The ``suffix`` determines the file
type. The following suffixes/file types are possible:

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

If you are submitting families then a .ped file for each family and an
all_samples ped file (cat *ped > all_samples.ped) must be supplied
within your working directory.  This is a tab-delimited file named
<family>.ped (where <family> is the family ID in the title of the
corresponding fastq files) with no header and one individual per line
according to the following pattern:

family_id sample_id father_id mother_id sex phenotype

family_id and sample_id should correspond to <family> and
<family>-<sample> in the sra/fastq filenames, father_id and mother_id
should be '0' if unknown, sex should be '1' if male, '2' if female and
'0' if unknown, and phenotype should be '1' if unaffected, '2' if
affected and '0' if unknown.

If you are running the functions to look for compound heterozygotes in
Multiplex families then there is a further requirement for the .ped
files.  The phasing tools expect a trio and therefore any other family
members (other than parents and one child) must be labelled as
unrelated.  That is, the first additional family member could be
labelled "family0" in the family_id column, and subsequent additional
family members could be "family1", "family2" and so on.  For example,
a multiplex family called Multiplex1 may have two parents and two
affected children.  The .ped file would look like this:

Multiplex1 ID1 0 0 1 1
Multiplex1 ID2 0 0 2 1
Multiplex1 ID3 ID1 ID2 1 2
Family0 ID4 ID1 ID2 2 2

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

Pipeline output
===============

The major output is a csvdb containing quality control information by
sample and variant information by family and an html report with
similar information.

Example
=======

ToDo: make exome sequencing example

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Requirements:

* BWA >= 0.7.8
* picardtools >= 1.106
* samtools >= 1.1
* GATK >= 2.7
* snpEff >= 4.0
* Gemini >= ?
* VCFtools >= 0.1.8a

Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *
from rpy2.robjects import r as R
import sys
import os
import csv
import glob
import re
import sqlite3
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.PipelineExome as PipelineExome
import PipelineExomeAncestry as PipelineExomeAncestry
import decimal
import pandas as pd
import numpy as np

###############################################################################
###############################################################################
###############################################################################
# load options from the config file


P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0], "pipeline.ini"])

PARAMS = P.PARAMS
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

matches = glob.glob("*.fastq.1.gz") + glob.glob(
    "*.fastq.gz") + glob.glob("*.sra") + glob.glob("*.csfasta.gz")


def getGATKOptions():
    return "-l mem_free=1.4G"

###############################################################################
###############################################################################
###############################################################################
# Load target and sample data into database
# The following functions are designed to upload meta-data to the csvdb
# These haven't been fully implemented yet


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["roi_bed"], "roi.load")
def loadROI(infile, outfile):
    '''Import regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    header = "chr,start,stop,feature"
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --header-names=%(header)s
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run()

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["roi_to_gene"], "roi2gene.load")
def loadROI2Gene(infile, outfile):
    '''Import genes mapping to regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run()

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@files(PARAMS["samples"], "samples.load")
def loadSamples(infile, outfile):
    '''Import sample information into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
            | cgat csv2db %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Alignment to a reference genome


@follows(mkdir("bam"))
@transform(INPUT_FORMATS, REGEX_FORMATS, r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA-MEM (output=SAM), convert to BAM,
    sort and index BAM file'''
    job_options = "-l mem_free=8G"
    job_threads = PARAMS["bwa_threads"]
    track = P.snip(os.path.basename(outfile), ".bam")
    m = PipelineMapping.BWAMEM(remove_unique=PARAMS["bwa_remove_non_unique"])
    statement = m.build((infiles,), outfile)
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Post-alignment QC


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.picard_stats")
def PicardAlignStats(infile, outfile):
    '''Run Picard CollectMultipleMetrics on each BAM file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    PipelineMappingQC.buildPicardAlignmentStats(infile, outfile, genome)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(PicardAlignStats, "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# GATK


@follows(loadPicardAlignStats, mkdir("gatk"))
@transform(mapReads, regex(r"bam/(\S+).bam"), r"gatk/\1.readgroups.bam")
def GATKReadGroups(infile, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = re.sub(r'-\w+-\w+\.bam', '', os.path.basename(infile))
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = PARAMS["gatk_threads"]
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    PipelineExome.GATKReadGroups(infile, outfile, genome,
                                 library, platform,
                                 platform_unit, track)
    IOTools.zapFile(infile)

###############################################################################
###############################################################################
###############################################################################
# Remove duplicates, realign and recalibrate lane-by-lane


@transform(GATKReadGroups,
           regex(r"gatk/(\S+).readgroups.bam"),
           r"gatk/\1.dedup.bam")
def RemoveDuplicatesLane(infile, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.buildPicardDuplicateStats(infile, outfile)
    IOTools.zapFile(infile)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(RemoveDuplicatesLane, "picard_duplicate_stats_lane.load")
def loadPicardDuplicateStatsLane(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)

###############################################################################


@transform(RemoveDuplicatesLane,
           regex(r"gatk/(\S+).dedup.bam"),
           r"gatk/\1.realigned2.bam")
def GATKIndelRealignLane(infile, outfile):
    '''realigns around indels using GATK'''
    threads = PARAMS["gatk_threads"]
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    PipelineExome.GATKIndelRealign(infile, outfile, genome, intervals, padding,
                                   threads)
    IOTools.zapFile(infile)

###############################################################################


@transform(GATKIndelRealignLane,
           regex(r"gatk/(\S+).realigned2.bam"),
           r"gatk/\1.bqsr.bam")
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    intrack = P.snip(os.path.basename(infile), ".bam")
    outtrack = P.snip(os.path.basename(outfile), ".bam")
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    if PARAMS["targetted"]:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(intrack + ".bai", outtrack + ".bai")
    else:
        PipelineExome.GATKBaseRecal(infile, outfile, genome, intervals,
                                    padding, dbsnp, solid_options)
    IOTools.zapFile(infile)

###############################################################################
###############################################################################
###############################################################################
# Merge BAMs across different lanes for the same sample


@collate(GATKBaseRecal,
         regex(r"gatk/(.*).bqsr.bam"),
         r"gatk/\1.merged.bam")
def mergeBAMs(infiles, outfile):
    '''merges BAMs for a single sample over multiple lanes'''
    inputfiles = " INPUT=".join(infiles)
    outf = open(outfile + ".count", "w")
    outf.write(str(len(infiles)))
    outf.close()
    statement = '''MergeSamFiles
                   INPUT=%(inputfiles)s
                   OUTPUT=%(outfile)s
                   ASSUME_SORTED=true; checkpoint; '''
    statement += '''samtools index %(outfile)s ;''' % locals()
    P.run()

    for inputfile in infiles:
        IOTools.zapFile(inputfile)

###############################################################################
###############################################################################
###############################################################################
# Remove duplicates sample-by-sample


@transform(mergeBAMs,
           regex(r"gatk/(\S+).merged.bam"),
           add_inputs(r"gatk/\1.merged.bam.count"),
           r"gatk/\1.dedup2.bam")
def RemoveDuplicatesSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    infile, countfile = infiles
    countf = open(countfile, "r")
    if countf.read() != '1':
        PipelineMappingQC.buildPicardDuplicateStats(infile, outfile)
    else:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(infile + ".bai", outfile + ".bai")
        IOTools.zapFile(infile)


###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(RemoveDuplicatesSample, "picard_duplicate_stats_sample.load")
def loadPicardDuplicateStatsSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# Realign sample-by-sample


@transform(RemoveDuplicatesSample,
           regex(r"gatk/(\S+).dedup2.bam"),
           add_inputs(r"gatk/\1.merged.bam.count"),
           r"gatk/\1.realigned.bam")
def GATKIndelRealignSample(infiles, outfile):
    '''realigns around indels using GATK'''
    infile, countfile = infiles
    threads = PARAMS["gatk_threads"]
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    countf = open(countfile, "r")
    if countf.read() > '1':
        PipelineExome.GATKIndelRealign(infiles[0], outfile, genome, intervals,
                                       padding, threads)
    else:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(infile + ".bai", outfile + ".bai")
#    IOTools.zapFile(infile)

###############################################################################
###############################################################################
###############################################################################
# Coverage of targetted area


@transform(GATKIndelRealignSample,
           regex(r"gatk/(\S+).realigned.bam"),
           r"gatk/\1.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed
    file using Picard'''
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    PipelineMappingQC.buildPicardCoverageStats(infile, outfile,
                                               baits, regions)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    '''Import coverage statistics into SQLite'''
    PipelineMappingQC.loadPicardCoverageStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# Guess sex


@follows(mkdir("xy_ratio"))
@transform(GATKIndelRealignSample,
           regex(r"gatk/(\S+).realigned.bam"),
           r"xy_ratio/\1.sex")
def calcXYratio(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads
    per megabase of sequence on X and Y'''
    PipelineExome.guessSex(infile, outfile)

###############################################################################


@merge(calcXYratio, "xy_ratio/xy_ratio.tsv")
def mergeXYRatio(infiles, outfile):
    '''merge XY ratios from all samples and load into database'''
    inlist = " ".join(infiles)
    statement = '''cgat combine_tables
                   --add-file-prefix --regex-filename="xy_ratio/(\S+).sex"
                   --no-titles --missing-value=0 --ignore-empty
                   -L %(outfile)s.log -v 6
                   --cat=Track %(inlist)s
                   > %(outfile)s'''
    P.run()

###############################################################################


@transform(mergeXYRatio, regex(r"xy_ratio/xy_ratio.tsv"),
           r"xy_ratio/xy_ratio.load")
def loadXYRatio(infile, outfile):
    '''load into database'''
    P.load(infile, outfile, "--header-names=Track,X,Y,XY_ratio")

###############################################################################
###############################################################################
###############################################################################
# Variant Calling


@follows(mkdir("variants"))
@transform(GATKIndelRealignSample, regex(r"gatk/(\S+).realigned.bam"),
           r"variants/\1.haplotypeCaller.g.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in individuals'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    job_options = getGATKOptions()
    job_threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    options = PARAMS["gatk_hc_options"]
    PipelineExome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options)

###############################################################################


@merge(haplotypeCaller, "variants/all_samples.vcf")
def genotypeGVCFs(infiles, outfile):
    '''Joint genotyping of all samples together'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    inputfiles = " --variant ".join(infiles)
    options = PARAMS["gatk_genotype_options"]
    PipelineExome.genotypeGVCFs(inputfiles, outfile, genome, options)

###############################################################################
###############################################################################
###############################################################################
# HapMap Genotypes


@follows(mkdir("hapmap"))
@files(PARAMS["hapmap_vcf"], "hapmap/hapmap_exome.bed")
def SelectExonicHapmapVariants(infile, outfile):
    '''Select variants from HapMap project that are located within
    exome target regions. Assumes bgzipped & tabix indexed Hapmap VCF file.'''
    bed = PARAMS["roi_regions"]
    statement = '''tabix -B %(infile)s %(bed)s |
                   awk '{OFS="\\t"; if (!/^#/){print $1,$2-1,$2}}'
                   > %(outfile)s''' % locals()
    P.run()

###############################################################################


@follows(SelectExonicHapmapVariants)
@transform(GATKIndelRealignSample,
           regex(r"gatk/(\S+).realigned.bam"),
           add_inputs("hapmap/hapmap_exome.bed"),
           r"hapmap/\1.hapmap.vcf")
def HapMapGenotype(infiles, outfile):
    '''Genotype HapMap SNPs using HaplotypeCaller in each individual'''
    infile, intervals = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    padding = PARAMS["hapmap_padding"]
    options = PARAMS["hapmap_hc_options"]
    PipelineExome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options)

###############################################################################


@transform(HapMapGenotype, regex(r"hapmap/(\S+).hapmap.vcf"),
           r"hapmap/\1.hapmap.vcf.gz")
def indexVCFs(infile, outfile):
    '''Genotype HapMap SNPs using HaplotypeCaller in each individual'''
    statement = '''bgzip -c %(infile)s > %(outfile)s;
                   tabix -p vcf %(outfile)s; '''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Compare Hapmap genotypes to assess relatedness


@follows(indexVCFs, mkdir("hapmap/vcfcompare"))
@permutations("hapmap/*.hapmap.vcf.gz",
              formatter("hapmap/(\S+).hapmap.vcf.gz$"),
              2,
              r"hapmap/vcfcompare/{basename[0][0]}_vs_{basename[1][0]}.vcfcompare")
def vcfCompare(infiles, outfile):
    '''Compare HapMap genotypes from each pair of calls '''
    sample1, sample2 = infiles
    name1 = P.snip(os.path.basename(sample1), ".hapmap.vcf.gz")
    name2 = P.snip(os.path.basename(sample2), ".hapmap.vcf.gz")
    statement = '''vcf-compare -g -m %(name1)s:%(name2)s
                   %(sample1)s %(sample2)s > %(outfile)s'''
    P.run()

###############################################################################


@transform(vcfCompare,
           regex(r"hapmap/vcfcompare/(\S+).vcfcompare"),
           r"hapmap/vcfcompare/\1.ndr")
def parseVcfCompare(infile, outfile):
    '''Extract non-reference discordance rate for each comparison'''
    statement = '''cat %(infile)s
                   | grep "Non-reference Discordance Rate (NDR):"
                   | cut -f 3
                   > %(outfile)s'''
    P.run()

###############################################################################


@merge(parseVcfCompare, "hapmap/vcfcompare/ndr.tsv")
def mergeNDR(infiles, outfile):
    '''Merge NDR values for all files'''
    outf = open(outfile, "w")
    outf.write("sample1\tsample2\tNDR\n")
    for infile in infiles:
        inbase = P.snip(os.path.basename(infile), ".ndr")
        sample1, sample2 = inbase.split("_vs_")
        inf = open(infile, "r")
        ndr = inf.readline().strip()
        result = "\t".join([sample1, sample2, ndr])
        outf.write(result + "\n")
        inf.close()
    outf.close()

###############################################################################


@transform(mergeNDR, regex(r"hapmap/vcfcompare/ndr.tsv"),
           r"hapmap/vcfcompare/ndr.load")
def loadNDR(infile, outfile):
    '''Load NDR into database'''
    P.load(infile, outfile)

###############################################################################
###############################################################################
###############################################################################
# Variant Annotation


@transform(genotypeGVCFs,
           regex(r"variants/all_samples.vcf"),
           r"variants/all_samples.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''snpEff.sh eff
                   -c %(config)s
                   -v %(snpeff_genome)s
                   -o gatk %(infile)s > %(outfile)s''' % locals()
    P.run()

###############################################################################


@transform(annotateVariantsSNPeff,
           regex(r"variants/all_samples.snpeff.vcf"),
           r"variants/all_samples.snpeff.table")
def vcfToTableSnpEff(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["annotation_snpeff_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(vcfToTableSnpEff, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadTableSnpEff(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")


###############################################################################
# GATK Variant Annotator


@merge(GATKIndelRealignSample, "gatk/all_samples.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for use in VQSR'''
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')

###############################################################################


@follows(annotateVariantsSNPeff)
@transform(genotypeGVCFs,
           regex(r"variants/all_samples.vcf"),
           add_inputs(listOfBAMs, r"variants/all_samples.snpeff.vcf"),
           r"variants/all_samples.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    vcffile, bamlist, snpeff_file = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    annotations = PARAMS["gatk_variant_annotations"]
    PipelineExome.variantAnnotator(vcffile, bamlist, outfile, genome,
                                   dbsnp, annotations, snpeff_file)


###############################################################################
# SNP Recalibration


@transform(variantAnnotator,
           regex(r"variants/all_samples.annotated.vcf"),
           r"variants/all_samples.snp_vqsr.recal")
def variantRecalibratorSnps(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    job_options = getGATKOptions()
    job_threads = PARAMS["gatk_threads"]
    track = P.snip(outfile, ".recal")
    kgenomes = PARAMS["gatk_kgenomes"]
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    mode = 'SNP'
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode, dbsnp,
                                      kgenomes, hapmap, omni)

###############################################################################


@follows(variantRecalibratorSnps)
@transform(variantAnnotator,
           regex(r"variants/all_samples.annotated.vcf"),
           add_inputs(r"variants/all_samples.snp_vqsr.recal",
                      r"variants/all_samples.snp_vqsr.tranches"),
           r"variants/all_samples.snp_vqsr.vcf")
def applyVariantRecalibrationSnps(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    mode = 'SNP'
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches,
                                            outfile, genome, mode)

###############################################################################
# Indel recalibration


@transform(applyVariantRecalibrationSnps,
           regex(r"variants/all_samples.snp_vqsr.vcf"),
           r"variants/all_samples.vqsr.recal")
def variantRecalibratorIndels(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    job_options = getGATKOptions()
    job_threads = PARAMS["gatk_threads"]
    track = P.snip(outfile, ".recal")
    mills = PARAMS["gatk_mills"]
    dbsnp = PARAMS["gatk_dbsnp"]
    mode = 'INDEL'
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode,
                                      mills=mills)

###############################################################################


@follows(variantRecalibratorIndels)
@transform(applyVariantRecalibrationSnps,
           regex(r"variants/all_samples.snp_vqsr.vcf"),
           add_inputs(r"variants/all_samples.vqsr.recal",
                      r"variants/all_samples.vqsr.tranches"),
           r"variants/all_samples.vqsr.vcf")
def applyVariantRecalibrationIndels(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    mode = 'INDEL'
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches,
                                            outfile, genome, mode)


###############################################################################
# SnpSift

@transform(applyVariantRecalibrationIndels,
           regex(r"variants/all_samples.vqsr.vcf"),
           r"variants/all_samples_dbnsfp.snpsift.vcf")
def annotateVariantsDBNSFP(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]
    dbNSFP = PARAMS["annotation_dbnsfp"]
    dbN_annotators = PARAMS["annotation_dbnsfpannotators"]
    if len(dbN_annotators) == 0:
        annostring = ""
    else:
        annostring = "-f %s" % dbN_annotators

    statement = """SnpSift.sh dbnsfp -db %(dbNSFP)s -v %(infile)s
                   %(annostring)s >
                   %(outfile)s;"""
    P.run()


@transform(annotateVariantsDBNSFP,
           regex(r"variants/all_samples_dbnsfp.snpsift.vcf"),
           r"variants/all_samples_clinvar.snpsift.vcf")
def annotateVariantsClinvar(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]
    clinvar = PARAMS["annotation_clinvar"]
    intout = outfile.replace("samples", "samples_clinvar")
    statement = """SnpSift.sh annotate %(clinvar)s
                %(infile)s > %(outfile)s;"""
    P.run()


@transform(annotateVariantsClinvar,
           regex(r"variants/all_samples_clinvar.snpsift.vcf"),
           r"variants/all_samples_exac.snpsift.vcf")
def annotateVariantsExAC(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]
    exac = PARAMS["annotation_exac"]
    statement = """SnpSift.sh annotate
                %(exac)s
                %(infile)s > %(outfile)s;"""
    P.run()


@transform(annotateVariantsExAC,
           regex(r"variants/all_samples_exac.snpsift.vcf"),
           r"variants/all_samples_gwasc.snpsift.vcf")
def annotateVariantsGWASC(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]

    gwas_catalog = PARAMS["annotation_gwas_catalog"]
    statement = """SnpSift.sh gwasCat -db %(gwas_catalog)s
                   %(infile)s > %(outfile)s;"""
    P.run()


@transform(annotateVariantsGWASC,
           regex(r"variants/all_samples_gwasc.snpsift.vcf"),
           r"variants/all_samples_phastcons.snpsift.vcf")
def annotateVariantsPhastcons(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]
    genomeind = "%s/%s.fa.fai" % (
        PARAMS['general_genome_dir'],
        PARAMS['general_genome'])
    phastcons = PARAMS["annotation_phastcons"]
    intout = outfile.replace("samples", "samples_phastc")
    statement = """ln -sf %(genomeind)s %(phastcons)s/genome.fai;
                   SnpSift.sh phastCons %(phastcons)s %(infile)s >
                   %(outfile)s;"""
    P.run()


@transform(annotateVariantsPhastcons,
           regex(r"variants/all_samples_phastcons.snpsift.vcf"),
           r"variants/all_samples_1000G.snpsift.vcf")
def annotateVariants1000G(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
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
        statement = """SnpSift.sh annotate
                       %(vcf)s
                       %(tempin)s > %(tempout)s;
                       mv %(tempout)s %(tempin)s"""
        P.run()

    shutil.move(tempin, outfile)


@transform(annotateVariants1000G,
           regex(r"variants/all_samples_1000G.snpsift.vcf"),
           r"variants/all_samples_dbsnp.snpsift.vcf")
def annotateVariantsDBSNP(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = PARAMS["annotation_threads"]

    dbsnp = PARAMS["annotation_dbsnp"]
    statement = """SnpSift.sh annotate
                %(dbsnp)s
                %(infile)s > %(outfile)s;"""

    P.run()


@follows(annotateVariantsDBSNP)
def annotateVariantsSNPsift():
    pass


@transform(annotateVariantsDBSNP,
           regex(r"variants/all_samples_dbsnp.snpsift.vcf"),
           r"variants/all_samples.vep.vcf")
def annotateVariantsVEP(infile, outfile):
    '''
    Adds annotations as specified in the pipeline.ini using Ensembl
    variant effect predictor (VEP).
    '''
    # infile - VCF
    # outfile - VCF with vep annotations
    job_options = "-l mem_free=6G"
    job_threads = 4

    VEP = PARAMS["annotation_vepannotators"].split(",")
    vep_annotators = PARAMS["annotation_vepannotators"]
    vep_path = PARAMS["annotation_veppath"]
    vep_cache = PARAMS["annotation_vepcache"]
    vep_species = PARAMS["annotation_vepspecies"]
    vep_assembly = PARAMS["annotation_vepassembly"]
    if len(vep_annotators) != 0:
        annostring = vep_annotators
        statement = '''perl %(vep_path)s/variant_effect_predictor.pl
                       --cache --dir %(vep_cache)s --vcf
                       --species %(vep_species)s
                       --fork 2
                       --assembly %(vep_assembly)s --input_file %(infile)s
                       --output_file %(outfile)s --force_overwrite
                       %(annostring)s --offline;'''
        P.run()


@follows(mkdir("variant_tables"))
@transform(GATKIndelRealignSample, regex(r"gatk/(.*).realigned.bam"),
           add_inputs(annotateVariantsVEP), r"variant_tables/\1.tsv")
def makeAnnotationsTables(infiles, outfile):
    '''
    Converts the multi sample vcf generated with Haplotype caller into
    a single table for each sample contain only positions called as variants
    in that sample and all columns from the vcf.
    '''
    bamname = infiles[0]
    inputvcf = infiles[1]
    TF = P.getTempFilename(".")
    samplename = bamname.replace(".realigned.bam",
                                 ".bam").replace("gatk/", "")

    statement = '''bcftools view -h %(inputvcf)s |
                   awk -F '=|,' '$1=="##INFO" || $1=="##FORMAT"
                   {printf("%%s\\t%%s\\n", $3, $9)}'
                   | sed 's/>//g' > %(TF)s'''
    P.run()
    cols = []
    colds = []
    cols2 = []
    for line in IOTools.openFile(TF).readlines():
        if line.split("\t")[0] != "Samples":
            cols.append("[%%%s]" % line.split("\t")[0])
            colds.append(line.split("\t")[1].strip().replace(" ", "_"))
            cols2.append(line.split("\t")[0])
    l1 = "\t".join(cols2)
    l2 = "\t".join(colds)
    out = open(outfile, "w")
    out.write('''CHROM\tPOS\tQUAL\tID\tFILTER\tREF1\tALT\tGT\t%s\n\
                 chromosome\tposition\tquality\tid\tfilter\tref\talt\t\
                 genotype\t%s\n''' % (l1, l2))
    out.close()
    cstring = "\\t".join(cols)
    cstring = "%CHROM\\t%POS\\t%QUAL\\t%ID\\t%FILTER\\t%REF\\t\
               %ALT\\t[%GT]\\t" + cstring
    if PARAMS['test'] == 1:
        statement = '''bcftools query
                   -f '%(cstring)s\\n'
                   -i 'FILTER=="PASS" && GT!="0/0" && GT!="./."'
                   %(inputvcf)s >> %(outfile)s'''
    else:
        statement = '''bcftools query -s %(samplename)s
                   -f '%(cstring)s\\n'
                   -i 'FILTER=="PASS" && GT!="0/0" && GT!="./."'
                   %(inputvcf)s >> %(outfile)s'''
    P.run()


###############################################################################
# Ancestry Functions #
HAPMAP = "%s/*txt.gz" % PARAMS['hapmap_loc']


@merge(HAPMAP, ["snpdict.json",
                "randomsnps.tsv"])
def makeRandomSNPSet(infiles, outfiles):
    '''
    Generates a random set of SNPs to use to characterise the ancestry
    and relatedness of the samples.

    1. Finds all SNPs which have a known genotype frequency for all 11
       hapmap ancestries
    2. Picks a random 50000 of these SNPs
    3. Stores a list of these SNPs as randomsnps.tsv
    4. Builds a dictionary of dictionaries of dictionaries where:
           At Level 1 each key is a HapMap ancestry ID
               each of these keys leads to another dictionary - level 2.
           At Level 2 each key is a dbSNP SNP ID
               each of these keys leads to another dictionary - level 3.
           At Level 3 each key is a genotype:
               the values of these keys are the genotype frequency of this
               genotype at this SNP in this ancestry.
       e.g. snpdict['ASW']['rs000000001']['CT'] -->
            frequency of the CT genotype at the rs0000001 SNP in the
            ASW population
    5.  Stores this dictionary in json format as randomsnps.json
    '''
    rs = PARAMS['general_randomseed']
    PipelineExomeAncestry.MakeSNPFreqDict(infiles, outfiles, rs, submit=True)


@follows(mkdir("sample_genotypes.dir"))
@transform(makeAnnotationsTables, regex("variant_tables/(.*).tsv"),
           add_inputs(makeRandomSNPSet), r"sample_genotypes.dir/\1.tsv")
def getSampleGenotypes(infiles, outfile):
    '''
    Fetches the genotype from the variant tables for all samples
    for SNPs in the hapmap sample from makeRandomSNPSet.

    Complex sites are ignored (as simple SNPs are sufficient for these
    calculations).
    These are:
        Sites which failed QC (column 3 in the variant table is not PASS)
        Sites with more than 2 alleles defined (column 6 in the variant table
        contains more than one alternative allele)
        SNPs with more than one ID
        Indels
    '''
    snps = set([line.strip()
                for line in IOTools.openFile(infiles[1][1]).readlines()])
    PipelineExomeAncestry.GenotypeSNPs(infiles[0], snps, outfile, submit=True)


@merge(getSampleGenotypes, "calledsnps.tsv")
def concatenateSNPs(infiles, outfile):
    '''
    Lists all the SNPs in the sample from makeRandomSNPSet where a variant has
    been called in any of the input files.
    Stores the reference genotype, chromosome and position to use later.
    '''
    snps = set()
    for f in infiles:
        with IOTools.openFile(f) as inp:
            for line in inp:
                line = line.strip().split("\t")
                snps.add("%s\t%s\t%s\t%s" % (line[0], line[4],
                                             line[1], line[2]))
    out = IOTools.openFile(outfile, "w")
    for snp in snps:
        out.write("%s\n" % (snp))
    out.close()


@follows(mkdir("sample_freqs"))
@transform(getSampleGenotypes, regex("sample_genotypes.dir/(.*).tsv"),
           add_inputs(concatenateSNPs, makeRandomSNPSet),
           [r"sample_freqs/\1.tsv", r"sample_freqs/\1_anc.tsv"])
def calculateAncestry(infiles, outfiles):
    '''
    Takes the data stored in MakeRandomSNPSet and the genotype of each sample
    at each site in calledsnps.tsv and tabulates the frequency of this
    genotype in each of the HapMap ancestry categories.
    The overall probability of each ancestry is then calculated as the
    product of these frequencies. These can only be used in comparison to
    each other - to show which of the 11 ancestries is most probable.
    '''
    calledsnps = infiles[1]
    snpdict = infiles[2][0]
    PipelineExomeAncestry.CalculateAncestry(infiles[0], calledsnps, snpdict,
                                            outfiles, submit=True)


@merge(calculateAncestry, "ancestry_estimate.tsv")
def mergeAncestry(infiles, outfile):
    '''
    Merges the output of the ancestry calculations for each sample into a
    single table.
    Draws a plot showing the score (probabilty) for each individual
    for their assigned ancestry and the second closest match.
    x = individual
    y = score
    large diamonds represent the best match and small triangles the second
    best match
    '''
    out = IOTools.openFile(outfile, "w")
    for f in infiles:
        infile = f[1]
        scores = []
        ancs = []
        for line in IOTools.openFile(infile).readlines():
            line = line.strip().split("\t")
            anc = line[0]
            score = decimal.Decimal(line[1])
            ancs.append(anc)
            scores.append(score)
        z = list(zip(ancs, scores))
        s = sorted(z, key=lambda x: x[1])[::-1]
        out.write("%s\t%s\t%s\t%s\t%s\n" % (f[0],
                                            s[0][0], s[0][1],
                                            s[1][0], s[1][1]))
    out.close()
    PipelineExomeAncestry.PlotAncestry(outfile)


@active_if(len(matches) > 1)
@merge((calculateAncestry, concatenateSNPs),
       ["all_samples.ped", "all_samples.map"])
def makePed(infiles, outfiles):
    '''
    Generates the required input for the Plink and King software packages.
       - PED file - columns are SNPs and rows are samples, each cell is the
         genotype of the sample at this SNP
       - MAP file - rows are SNPs in the same order as the columns in the ped
         file, each row shows chromosome, snp id and position.
    '''
    PipelineExomeAncestry.MakePEDFile(infiles, outfiles)


@active_if(len(matches) > 1)
@transform(makePed, suffix(".ped"), ".ibs0")
def runPlinkandKing(infiles, outfile):
    '''
    Uses King software to calculate relatedness of pairs of samples as
    described here - http://people.virginia.edu/~wc9c/KING/manual.html
    Plink is used just to reformat the PED files into the format required by
    King.
    '''
    pref = infiles[0].split(".")[0]
    k = PARAMS['king_path']
    p = PARAMS['king_plink']
    statement = """
    %(p)s/plink --file %(pref)s --make-bed --no-fid --noweb --map3
          --no-parents --no-sex --no-pheno;
    %(k)s/king -b plink.bed --binary --prefix %(pref)s;
    %(k)s/king -b %(pref)s.bgeno --kinship --ibs --prefix %(pref)s"""
    P.run()


@active_if(len(matches) > 1)
@merge(runPlinkandKing, "family_estimate.tsv")
def calculateFamily(infile, outfile):
    '''
    Translates and filters the output from King.
    Pairs of related samples are written to the output file along with the
    degree of relatedness.  Degrees are decided using thresholds from
    the King documentation, here
    http://people.virginia.edu/~wc9c/KING/manual.html
    '''
    PipelineExomeAncestry.CalculateFamily(infile, outfile)


@active_if(len(matches) > 1)
@merge(makeAnnotationsTables, "sex_estimate.tsv")
def calculateSex(infiles, outfile):
    '''
    Approximates the sex of the sample using the data in the variant table.
    Basic estimate based on heterozygosity on the X chromosome - genotypes
    are taken for all called variants on X passing QC and the percentage
    of heterozygotes is taken.
    This tends to produce two clear populations so Kmeans clustering
    is used to split the data into two - male and female.  Samples which are
    unclear are marked in the output.
    '''
    PipelineExomeAncestry.CalculateSex(infiles, outfile, submit=True)


@merge((mergeAncestry, calculateFamily, calculateSex), "summary.tsv")
def summarise(infiles, outfile):
    '''
    Builds a summary table for ancestry, family and sex.
    '''
    ancestry = pd.read_csv(infiles[0], sep="\t", header=None)
    if len(matches) > 1:
        family = pd.read_csv(infiles[1], sep="\t", header=None)
        sex = pd.read_csv(infiles[2], sep="\t", header=None)
        ancestry[0] = [a[-1] for a in ancestry[0].str.split("/")]
        family[0] = [a[-1] for a in family[0].str.split("/")]
        sex[0] = [a[-1] for a in sex[0].str.split("/")]
        sex = sex.drop([1, 2], 1)
        sex.columns = ['id', 'sex', 'sex_significance']
        family.columns = ['id', 'related_to', 'degree_relatedness']
        family['degree_relatedness'] = family['degree_relatedness'].astype(int)
        ancestry = ancestry.drop([2, 3, 4], 1)
        ancestry.columns = ['id', 'ancestry']
        summary = ancestry.merge(sex).merge(family, 'left')
        summary = summary.fillna("NA")
        summary.to_csv(outfile, sep="\t")
    else:
        ancestry[0] = [a[-1] for a in ancestry[0].str.split("/")]
        ancestry = ancestry.drop([2, 3, 4], 1)
        ancestry.columns = ['id', 'ancestry']
        ancestry = ancestry.fillna("NA")
        ancestry.to_csv(outfile, sep="\t")


@follows(summarise)
def ancestry():
    pass


###############################################################################
# Quality Filtering #

@follows(mkdir("variant_tables_highqual"))
@transform(makeAnnotationsTables, regex("variant_tables/(.*).tsv"),
           [r'variant_tables_highqual/\1.tsv',
            r'variant_tables_highqual/\1_failed.tsv'])
def qualityFilterVariants(infile, outfiles):
    '''
    Filter variants based on quality.  Columns to filter on and
    how they should be filtered can be specified in the pipeline.ini.
    Currently only implemented to filter numeric columns.  "." is assumed
    to mean pass.
    '''
    qualstring = PARAMS['filtering_quality']
    qualfilter = PARAMS['filtering_quality_ft']
    PipelineExome.filterQuality(infile, qualstring, qualfilter,
                                outfiles, submit=True)


@follows(mkdir("variant_tables_rare"))
@transform(qualityFilterVariants, regex("variant_tables_highqual/(.*).tsv"),
           [r'variant_tables_rare/\1.tsv',
            r'variant_tables_rare/\1_failed.tsv'])
def rarityFilterVariants(infiles, outfiles):
    '''
    Filter out variants which are common in any of the exac or other
    population datasets as specified in the pipeline.ini.
    '''
    infile = infiles[0]
    thresh = PARAMS['filtering_rarethresh']
    freqs = PARAMS['filtering_freqs']
    exac = PARAMS['filtering_exac']
    PipelineExome.filterRarity(infile, exac, freqs, thresh, outfiles,
                               submit=True)


@follows(mkdir("variant_tables_damaging"))
@transform(rarityFilterVariants, regex("variant_tables_rare/(.*).tsv"),
           [r'variant_tables_damaging/\1.tsv',
            r'variant_tables_damaging/\1_failed.tsv'])
def damageFilterVariants(infiles, outfiles):
    '''
    Filter variants which have not been assessed as damaging by any
    of the specified tools.
    Tools and thresholds can be specified in the pipeline.ini.

    Does not account for multiple alt alleles - if any ALT allele has
    been assessed as damaging with any tool the variant is kept,
    regardless of if this is the allele called in the sample.

    '''
    infile = infiles[0]
    PipelineExome.filterDamage(infile, PARAMS['filtering_damage'], outfiles,
                               submit=True)


@follows(mkdir("variant_tables_family"))
@transform(damageFilterVariants, regex("variant_tables_damaging/(.*).tsv"),
           add_inputs(calculateFamily),
           [r'variant_tables_family/\1.tsv',
            r'variant_tables_family/\1_failed.tsv'])
def familyFilterVariants(infiles, outfiles):
    '''
    Filter variants according to the output of calculateFamily -
    only variants shared by both members of a family will be kept.
    '''
    if len(matches) > 1:
        infile = infiles[0][0]

        infilenam = infile.split("/")[-1]
        infilestem = "/".join(infile.split("/")[:-1])

        # figure out who is related to who
        families = [line.strip().split("\t")[:2]
                    for line in IOTools.openFile(infiles[1]).readlines()]
        infam = [line[0] for line in families] + [line[1] for line in families]

        # no relatives - copy the input file to the output file and generate
        # a blank "failed" file
        if infilenam not in infam or PARAMS['filtering_family'] == 0:
            shutil.copy(infile, outfiles[0])
            o = IOTools.openFile(outfiles[1], "w")
            o.close()
        else:
            for line in families:
                if infilenam in line:
                    i = line.index(infilenam)
                    if i == 0:
                        infile2 = "%s/%s" % (infilestem, line[1])
                    else:
                        infile2 = "%s/%s" % (infilestem, line[0])
            PipelineExome.filterFamily(infile, infile2, outfiles)
    else:
        infile = infiles[0][0]
        shutil.copy(infile, outfiles[0])
        out = IOTools.openFile(outfiles[1], "w")
        out.close()


@follows(familyFilterVariants)
def filter():
    pass


@follows(mkdir("genelists"))
@transform(familyFilterVariants, regex("variant_tables_family/(.*).tsv"),
           [r'genelists/\1.tsv',
            r'genelists/\1.bed'])
def makeGeneLists(infiles, outfiles):
    infile = infiles[0]
    outfile = outfiles[1]
    genes = PARAMS['general_geneset']

    # four %s because they need to be escaped in generating the statement
    # then again when submitting the P.run()
    statement = '''awk 'NR > 2 {printf("%%%%s\\t%%%%s\\t%%%%s\\n",\
    $1, $2, $2 + 1)}'\
    %(infile)s |\
    bedtools intersect -wo -a stdin -b %(genes)s > %(outfile)s''' % locals()
    P.run()

    geneids = set()
    with IOTools.openFile(outfile) as inp:
        for line in inp:
            line = line.strip().split("\t")
            details = line[11].split(";")
            for detail in details:
                r = re.search('gene_id', detail)
                if r:
                    geneid = detail.split(" ")[-1]
                    geneids.add(geneid.replace("\"", ""))
    out = IOTools.openFile(outfiles[0], "w")
    for geneid in geneids:
        out.write("%s\n" % geneid)
    out.close()


@follows(mkdir('final_variant_tables'))
@collate((makeGeneLists, familyFilterVariants), regex('(.*)/(.*).tsv'),
         r'final_variant_tables/\2.tsv')
def finalVariantTables(infiles, outfile):
    genes = infiles[0]
    variants = infiles[1][0]
    cols = PARAMS['filtering_columns'].split(",")
    PipelineExome.CleanVariantTables(genes, variants, cols, outfile,
                                     submit=True)

###############################################################################
###############################################################################
###############################################################################
# Genes of interest


@active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
@transform((annotateVariantsVEP),
           regex(r"variants/all_samples.snpsift.vcf"),
           r"variants/all_samples.genes.vcf")
def findGenes(infile, outfile):
    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
    if variant is within a gene of interest as defined in the ini
    file'''
    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
    statement = '''GenomeAnalysisTK -T VariantFiltration
                   -R %%(genome_dir)s/%%(gatkgenome)s.fa
                   --variant %(infile)s
                   --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
                   --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Tabulation


TABULATION_INPUT = {0: annotateVariantsVEP, 1: findGenes}


@transform(TABULATION_INPUT[PARAMS["annotation_add_genes_of_interest"]],
           regex(r"variants/all_samples.(snpsift|genes).vcf"),
           r"variants/all_samples.snpsift.table")
def vcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(vcfToTable, regex(r"variants/all_samples.snpsift.table"),
           r"variants/all_samples.snpsift.table.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")

##############################################################################
###############################################################################
###############################################################################
# Confirm parentage (do novo trios only)


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.parentage")
def confirmParentage(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
        'family', 'sample', 'father', 'mother', 'sex', 'status'])
    trio = P.snip(os.path.basename(pedfile), ".ped")
    trio = trio.replace(".", "_").replace("-", "_")
    database = PARAMS["database_name"]
    proband = None
    mother = None
    father = None
    for row in pedigree:
        if row['status'] == '2':
            proband = row['sample'].replace(".", "_").replace("-", "_")
            mother = row['mother'].replace(".", "_").replace("-", "_")
            father = row['father'].replace(".", "_").replace("-", "_")
    E.info("proband:" + proband)
    E.info("mother:" + mother)
    E.info("father:" + father)

    query = '''SELECT '%(trio)s' as family,
               hom_nonref_trio/(hom_nonref_parents+0.0) as hom_nonref_conc,
               child_mat_het_nonref/(maternal_hom_nonref+0.0) as maternal_conc,
               child_pat_het_nonref/(paternal_hom_nonref+0.0) as paternal_conc,
                *
               FROM
               (SELECT count(*) as hom_nonref_trio
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as hom_nonref_parents
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as maternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_mat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as paternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_pat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10)''' % locals()
    statement = '''sqlite3 %(database)s "%(query)s" > %(outfile)s
                   2> %(outfile)s.log''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# De novos


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.filtered.vcf")
def deNovoVariants(infiles, outfile):
    '''Filter de novo variants based on provided jexl expression'''
    job_options = getGATKOptions()
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(father)s").getDP()>=10&&vc.getGenotype("%(mother)s").getDP()>=10&&vc.getGenotype("%(child)s").getPL().0>20&&vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(child)s").getPL().2>0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(father)s").getPL().1>20&&vc.getGenotype("%(father)s").getPL().2>20&&vc.getGenotype("%(mother)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().1>20&&vc.getGenotype("%(mother)s").getPL().2>20&&vc.getGenotype("%(child)s").getAD().1>=3&&((vc.getGenotype("%(child)s").getAD().1)/(vc.getGenotype("%(child)s").getDP().floatValue()))>=0.25&&(vc.getGenotype("%(father)s").getAD().1==0||(vc.getGenotype("%(father)s").getAD().1>0&&((vc.getGenotype("%(father)s").getAD().1)/(vc.getGenotype("%(father)s").getDP().floatValue()))<0.05))&&(vc.getGenotype("%(mother)s").getAD().1==0||(vc.getGenotype("%(mother)s").getAD().1>0&&((vc.getGenotype("%(mother)s").getAD().1)/(vc.getGenotype("%(mother)s").getDP().floatValue()))<0.05))&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    PipelineExome.selectVariants(infile, outfile, genome, select)


@transform(deNovoVariants,
           regex(r"variants/(\S+).filtered.vcf"),
           r"variants/\1.filtered.table")
def tabulateDeNovos(infile, outfile):
    '''Tabulate de novo variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateDeNovos,
           regex(r"variants/(\S+).filtered.table"),
           r"variants/\1.filtered.table.load")
def loadDeNovos(infile, outfile):
    '''Load de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################


@follows(loadVariantAnnotation)
@transform("*Trio*.ped",
           regex(r"(\S*Trio\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.denovos.vcf")
def lowerStringencyDeNovos(infiles, outfile):
    '''Filter lower stringency de novo variants based on provided jexl
    expression'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().0==0&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    PipelineExome.selectVariants(infile, outfile, genome, select)


@transform(lowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.vcf"),
           r"variants/\1.denovos.table")
def tabulateLowerStringencyDeNovos(infile, outfile):
    '''Tabulate lower stringency de novo variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateLowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.table"),
           r"variants/\1.denovos.table.load")
def loadLowerStringencyDeNovos(infile, outfile):
    '''Load lower stringency de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Dominant


@follows(loadVariantAnnotation)
@transform("*Multiplex*.ped",
           regex(r"(\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.dominant.vcf")
def dominantVariants(infiles, outfile):
    '''Filter variants according to autosomal dominant disease model'''
    pedfile, infile = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
        if row['status'] == '1':
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
    affecteds_exp2 = '").getAD().1>=1&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
            '").isHomRef()'
    # for some weird reason the 1000G filter doesn't work on these particular
    # files - will add later when I've figured out what's wrong
    # currently 1000G filter is performed at the report stage (not in csvdb)
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().1==0&&vc.getGenotype("%(affecteds_exp2)s").getAD().1>=1%(unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    PipelineExome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(dominantVariants,
           regex(r"variants/(\S+).dominant.vcf"),
           r"variants/\1.dominant.table")
def tabulateDoms(infile, outfile):
    '''Tabulate dominant disease candidate variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateDoms, regex(r"variants/(\S+).dominant.table"),
           r"variants/\1.dominant.table.load")
def loadDoms(infile, outfile):
    '''Load dominant disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Recessive


@follows(loadVariantAnnotation)
@transform("*.ped",
           regex(
               r"(\S*Trio\S+|\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.recessive.vcf")
def recessiveVariants(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    parents = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
            if row['father'] != '0':
                parents += [row['father']]
            if row['mother'] != '0':
                parents += [row['mother']]
        elif row['status'] == '1' and row['sample'] not in parents:
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)) + \
            '").getPL().2!=0'
    if len(parents) == 0:
        parents_exp = ''
    else:
        parents_exp = '&&vc.getGenotype("' + \
            ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
            '").getPL().1==0'
    # need a way of specifying that other unaffecteds eg. sibs can't be
    # homozygous for alt allele
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().2==0%(unaffecteds_exp)s%(parents_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    PipelineExome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(recessiveVariants,
           regex(r"variants/(\S+).recessive.vcf"),
           r"variants/\1.recessive.table")
def tabulateRecs(infile, outfile):
    '''Tabulate potential homozygous recessive disease variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateRecs, regex(r"variants/(\S+).recessive.table"),
           r"variants/\1.recessive.table.load")
def loadRecs(infile, outfile):
    '''Load homozygous recessive disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# X-linked


@follows(loadVariantAnnotation)
@transform("*.ped", regex(r"(\S*Trio\S+|\S*Multiplex\S+).ped"),
           add_inputs(r"variants/all_samples.snpsift.vcf"),
           r"variants/\1.xlinked.vcf")
def xlinkedVariants(infiles, outfile):
    '''Find maternally inherited X chromosome variants in male patients'''
    track = P.snip(os.path.basename(outfile), ".vcf")
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    pedfile, infile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t',
                              fieldnames=['family', 'sample', 'father',
                                          'mother', 'sex', 'status'])
    affecteds = []
    mothers = []
    male_unaffecteds = []
    female_unaffecteds = []
    for row in pedigree:
        if row['status'] == '2' and row['mother'] != '0':
            if row['mother'] not in mothers:
                mothers += [row['mother']]
            affecteds += [row['sample']]
        elif row['status'] == '1' and row['sex'] == '1':
            male_unaffecteds += [row['sample']]
        elif row['status'] == '1' and row['sex'] != '1' and row['sample'] not in mothers:
            female_unaffecteds += [row['sample']]
    if len(affecteds) == 0:
        affecteds_exp = ''
    else:
        affecteds_exp = '").isHomRef()&&!vc.getGenotype("'.join(affecteds)
#        affecteds_exp = '&&!vc.getGenotype("' + \
#        ('").isHomRef()&&!vc.getGenotype("'.join(affecteds)) + \
#        '").isHomRef()'
    if len(mothers) == 0:
        mothers_exp = ''
    else:
        mothers_exp = '&&vc.getGenotype("' + \
                      ('").getPL().1==0&&vc.getGenotype("'.join(mothers)) + \
                      '").getPL().1==0'
    if len(male_unaffecteds) == 0:
        male_unaffecteds_exp = ''
    else:
        male_unaffecteds_exp = '&&vc.getGenotype("' + \
                               ('").isHomRef()&&vc.getGenotype("'.join(
                                   male_unaffecteds)) + \
                               '").isHomRef()'
    if len(female_unaffecteds) == 0:
        female_unaffecteds_exp = ''
    else:
        female_unaffecteds_exp = '&&vc.getGenotype("' + \
                                 ('").getPL().2!=0&&vc.getGenotype("'.join(
                                     female_unaffecteds)) + \
            '").getPL().2!=0'
    select = '''CHROM=="X"&&!vc.getGenotype("%(affecteds_exp)s").isHomRef()%(mothers_exp)s%(male_unaffecteds_exp)s%(female_unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    PipelineExome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(xlinkedVariants,
           regex(r"variants/(\S+).xlinked.vcf"),
           r"variants/\1.xlinked.table")
def tabulateXs(infile, outfile):
    '''Tabulate potential X-linked disease variants'''
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(tabulateXs, regex(r"variants/(\S+).xlinked.table"),
           r"variants/\1.xlinked.table.load")
def loadXs(infile, outfile):
    '''Load X-linked disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Compound hets


# why does this not work from snpsift VCF file? DS
# @transform(annotateVariantsSNPsift,
#    regex(r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
@transform(annotateVariantsSNPeff,
           regex(
               r"variants/all_samples.snpeff.vcf"),
           add_inputs(r"all_samples.ped", r"gatk/all_samples.list"),
           r"variants/all_samples.phased.vcf")
def phasing(infiles, outfile):
    '''phase variants with GATK'''
    infile, pedfile, bamlist = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    statement = '''GenomeAnalysisTK -T PhaseByTransmission
                   -R %(genome)s
                   -V %(infile)s
                   -ped %(pedfile)s
                   -mvf %(infile)s.mvf
                   -o %(outfile)s ;''' % locals()
    P.run()

###############################################################################


@transform(phasing, regex(r"variants/all_samples.phased.vcf"),
           add_inputs(r"gatk/all_samples.list"),
           r"variants/all_samples.rbp.vcf")
def readbackedphasing(infiles, outfile):
    '''phase variants with ReadBackedPhasing'''
    job_memory = "32G"
    infile, bamlist = infiles
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    statement = '''GenomeAnalysisTK -T ReadBackedPhasing -nt 4
                   -R %(genome)s
                   -I %(bamlist)s
                   -V %(infile)s
                   -o %(outfile)s; ''' % locals()
    P.run()

###############################################################################


@follows(readbackedphasing)
@transform("*.ped",
           regex(r"(\S*Multiplex\S+|\S*Trio\S+).ped"),
           add_inputs(r"no_multiallelic_all_samples.rbp.vcf",
                      r"all_samples.ped"),
           r"variants/\1.compound_hets.table")
def compoundHets(infiles, outfile):
    '''Identify potentially pathogenic compound heterozygous variants
    (phasing with GATK followed by compound het calling using Gemini'''
    family, infile, pedfile = infiles
    family_id = P.snip(os.path.basename(family), ".ped")
    statement = '''gemini load -v %(infile)s
                   -p %(pedfile)s -t snpEff %(family_id)s.db ;
                   gemini comp_hets
                   --families %(family_id)s
                   --columns "chrom, start, end, ref, alt, codon_change, gene, qual, depth"
                   --filter
                   "(impact_severity = 'HIGH' OR impact_severity = 'MED')
                   AND (in_esp = 0 OR aaf_esp_all < 0.01)
                   AND (in_1kg = 0 OR aaf_1kg_all < 0.01)"
                   %(family_id)s.db > %(outfile)s;'''
    # rm -f %(family_id)s.db'''
    P.run()

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform(compoundHets, regex(r"variants/(\S+).compound_hets.table"),
           r"variants/\1.compound_hets.table.load")
def loadCompoundHets(infile, outfile):
    '''Load compound heterozygous variants into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# coverage over candidate genes


@active_if(PARAMS["coverage_calculate"] == 1)
@transform(listOfBAMs, regex(r"gatk/all_samples.list"), r"candidate.sample_interval_summary")
def candidateCoverage(infile, outfile):
    '''Calculate coverage over exons of candidate genes'''
    all_exons = PARAMS["coverage_all_exons"]
    candidates = PARAMS["coverage_candidates"]
    candidates = candidates.replace(",", " -e ")
    genome = PARAMS["genome_dir"] + "/" + PARAMS["gatkgenome"] + ".fa"
    threshold = PARAMS["coverage_threshold"]
    statement = '''zcat %(all_exons)s | grep -e %(candidates)s
                   | awk '{print $1 ":" $4 "-" $5}' - | sed 's/chr//' - >
                   candidate.interval_list ; ''' % locals()
    statement += '''zcat %(all_exons)s | grep -e %(candidates)s
                    | awk '{print $16}' - | sed 's/"//g;s/;//g' - >
                    candidate_gene_names.txt ;''' % locals()
    statement += '''GenomeAnalysisTK -T DepthOfCoverage -R %(genome)s
                    -o candidate -I %(infile)s -ct %(threshold)s -L candidate.interval_list ;''' % locals()
    P.run()

###############################################################################


@active_if(PARAMS["coverage_calculate"] == 1)
@transform(candidateCoverage, regex(r"candidate.sample_interval_summary"), r"candidate_coverage_plot.pdf")
def candidateCoveragePlots(infile, outfile):
    '''Produce plots of coverage'''
    rscript = PARAMS["coverage_rscript"]
    threshold = PARAMS["coverage_threshold"]
    statement = '''Rscript %(rscript)s %(infile)s candidate_gene_names.txt %(threshold)s %(outfile)s ;'''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# vcf statistics


@transform((annotateVariantsVEP), regex(
    r"variants/all_samples.snpsift.vcf"), r"variants/all_samples.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    statement = '''vcf-stats %(infile)s > %(outfile)s
                   2>>%(outfile)s.log;''' % locals()
    P.run()

###############################################################################


@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    E.info("Loading vcf stats...")
    statement = '''cgat vcfstats2db %(filenames)s >>
                   %(outfile)s; '''
    statement += '''cat vcfstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=vcf_stats >> %(outfile)s; '''
    statement += '''cat sharedstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=vcf_shared_stats >> %(outfile)s; '''
    statement += '''cat indelstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=indel_stats >> %(outfile)s; '''
    statement += '''cat snpstats.txt | cgat csv2db
                    %(csv2db_options)s --allow-empty-file --add-index=track
                    --table=snp_stats >> %(outfile)s; '''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Targets


@follows(loadVariantAnnotation,
         finalVariantTables,
         ancestry)
def testFromVariantRecal():
    pass


@follows(loadROI,
         loadROI2Gene,
         loadSamples)
def loadMetadata():
    pass


@follows(mapReads,
         loadPicardAlignStats)
def mapping():
    pass


@follows(GATKBaseRecal,
         loadPicardDuplicateStatsLane,
         GATKIndelRealignSample,
         loadCoverageStats)
def gatk():
    pass


@follows(loadXYRatio,
         indexVCFs,
         loadNDR)
def sampleFeatures():
    pass


@follows(haplotypeCaller,
         genotypeGVCFs)
def callVariants():
    pass


@follows(loadTableSnpEff,
         listOfBAMs,
         loadVariantAnnotation,
         finalVariantTables)
def annotation():
    pass


@follows(loadDeNovos)
def denovo():
    pass


@follows(loadLowerStringencyDeNovos)
def denovo2():
    pass


@follows(loadDoms)
def dominant():
    pass


@follows(loadRecs)
def recessive():
    pass


@follows(loadXs)
def xlinked():
    pass


@follows(loadCompoundHets)
def compoundHet():
    pass


@follows(denovo,
         dominant,
         recessive,
         compoundHet)
def filtering():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(mapping,
         gatk,
         # sampleFeatures,
         callVariants,
         annotation,
         filtering,
         ancestry,
         makeAnnotationsTables)
def full():
    pass

###############################################################################
###############################################################################
###############################################################################
# Reports


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
