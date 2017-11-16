"""===========================
Ancestral repeats pipeline
===========================

:Author: Andreas Heger
:Release: $Id: pipeline_ancestral_repeats.py 2876 2010-03-27 17:42:11Z andreas $
:Date: |today|
:Tags: Python

The ancestral repeats pipeline defines ancestral repeats for a pair of genomes
and computes rates for these.

This pipeline performs the following actions:

   * collect repeatmasker annotation from external
     databases. Currently implemented are:
      * UCSC
      * Ensembl
   * build pairwise genomic alignment from axt or maf files
   * define ancestral repeats
   * compute rates of ancestral repeats

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline expects a :term:`query` and :term:`target` genome. These
should be set in the general section.  For each genome there should
then be section on how to obtain the repeatmasker tracks. The default
configuration file gives an example.

Input
-----

The pipeline starts from an empty working directory. It will collect
the input data from directories specified in the configuration files.

The genomic alignment can both be build from :term:`axt` formatted
pairwise alignments and from :term:`maf` formatted multiple
alignments.

:term:`axt` formatted files (such as file:`chr1.hg19.mm10.net.axt.gz`)
are build by hierarchically from chains by selecting the
highest-scoring non-overlapping chains on top and then filling in the
gaps with lower scoring chains. A net is single-coverage for target
but not for query, unless it has been filtered to be single-coverage
on both target and query. Because it's single-coverage in the target,
it's no longer symmetrical. For ancestral repeat determination, use
:term:`axt` files that have been filtered to be single-coverage on
both query and target. By convention, the UCSC adds "rbest" to the net
filename in that case, such as: `hg19.panTro3.rbest.net.gz`.

:term:`maf` files currently only work if the :term:`query` genome is
the reference species in the maf files. This is a consequence of
:file:`maf2Axt` requiring that the strand of the reference species is
always positive and I have not figured out how to invert maf
alignments.

.. note::
   ENSEMBL import is not thoroughly tested.
   :term:`maf` formatted import is not thoroughly tested.

Type::

   python pipeline_ancestral_repeats.py --help

for command line help.

Requirements
------------



Output
======

The pipeline builds the following files:

aligned_repeats.psl.gz
   :term:`psl` formatted files of alignments between ancestral repeats

aligned_repeats.rates.gz
   rates between ancestral repeats

alignment.psl.gz
   :term:`psl` formatted genomic alignment between query and target.

<query>_rates.gff.gz
   :term:`gff` formatted file of ancestral repeats on the query. The score field is set
   to the estimated substitution rate of the repeat.

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_ancestral_repeats.tgz.

To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_ancestral_repeats.tgz
   tar -xvzf pipeline_ancestral_repeats.tgz
   cd pipeline_ancestral_repeats
   python <srcdir>/pipeline_ancestral_repeats.py make full

The example data builds ancestral repeats between human hg19:chr19 and
mouse mm9:chr7.

Code
====

"""
import sys
import os
import CGAT.Experiment as E
import logging as L
import CGATPipelines.PipelineUCSC as PipelineUCSC

from ruffus import *

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'query': "",
        'target': ""})
PARAMS = P.PARAMS

if os.path.exists("pipeline_conf.py"):
    L.info("reading additional configuration from pipeline_conf.py")
    exec(compile(open("pipeline_conf.py").read(), "pipeline_conf.py", 'exec'))


def getGenomes():
    '''return genome names of query and target.'''

    genome_query = os.path.join(PARAMS["genome_dir"], PARAMS["query"])
    genome_target = os.path.join(PARAMS["genome_dir"], PARAMS["target"])
    return genome_query, genome_target


@files([("%s/%s.idx" % (PARAMS["genome_dir"], x), "%s.sizes" % x)
        for x in (PARAMS["query"], PARAMS["target"])])
def buildSizes(infile, outfile):
    '''extract size information from genomes.'''
    outf = open(outfile, "w")
    for line in open(infile):
        data = line[:-1].split("\t")
        if len(data) >= 4:
            contig = data[0]
            outf.write("%s\t%s\n" % (contig, data[3]))
    outf.close()


if "axt_dir" in PARAMS:
    # build pairwise alignment from axt formatted data.'''
    @follows(buildSizes)
    @merge("%s/*.axt.gz" % PARAMS["axt_dir"],
           PARAMS["interface_alignment_psl"])
    def buildGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from axt files.'''

        try:
            os.remove(outfile)
        except OSError:
            pass

        for infile in infiles:
            E.info("adding %s" % infile)
            statement = '''gunzip < %(infile)s
            | axtToPsl
            /dev/stdin
            %(query)s.sizes
            %(target)s.sizes
            /dev/stdout
            | pslSwap /dev/stdin /dev/stdout
            | gzip >> %(outfile)s
            '''
            P.run()


elif "maf_dir" in PARAMS:
    @follows(buildSizes)
    @merge("%s/*.maf.gz" % PARAMS["maf_dir"], "alignment.raw.psl.gz")
    def buildRawGenomeAlignment(infiles, outfile):
        '''build pairwise genomic aligment from maf files.
        '''

        try:
            os.remove(outfile)
        except OSError:
            pass

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile:
                continue

            E.info("adding %s" % infile)

            genome_query, genome_target = getGenomes()

            statement = '''gunzip < %(infile)s
             | cgat maf2psl
                  --query=%(maf_name_query)s
                  --target=%(maf_name_target)s
                  --log=%(outfile)s.log
             | cgat psl2psl
                  --method=filter-fasta
                  --method=sanitize
                  --queries-tsv-file=%(genome_query)s
                  --target-psl-file=%(genome_target)s
                  --log=%(outfile)s.log
             | gzip
             >> %(outfile)s
             '''
            P.run()

    @transform(buildRawGenomeAlignment,
               suffix(".raw.psl.gz"),
               ".psl.gz")
    def buildGenomeAlignment(infile, outfile):
        '''remove non-unique alignments in genomic infile.'''

        statement = '''gunzip < %(infile)s
        | sort -k10,10 -k12,12n
        | cgat psl2psl
        --method=remove-overlapping-query
        --log=%(outfile)s.log
        | sort -k14,14 -k16,16n
        | cgat psl2psl
        --method=remove-overlapping-target
        --log=%(outfile)s.log
        | gzip
        >> %(outfile)s
        '''
        P.run()

    @follows(buildSizes)
    @merge("%s/*.maf.gz" % PARAMS["maf_dir"], PARAMS["interface_alignment_psl"])
    def buildGenomeAlignmentUCSCTools(infiles, outfile):
        '''build pairwise genomic aligment from maf files.'''

        try:
            os.remove(outfile)
        except OSError:
            pass

        for infile in infiles:
            # skip maf files without Hsap on top.
            if "other" in infile or "supercontig" in infile:
                continue

            E.info("adding %s" % infile)

            genome_query, genome_target = getGenomes()

            statement = '''gunzip < %(infile)s
            | mafToAxt
                  /dev/stdin
                  %(maf_name_target)s
                  %(maf_name_query)s
                  /dev/stdout
                  -stripDb
             | axtToPsl
                  /dev/stdin
                  %(target)s.sizes
                  %(query)s.sizes
                  /dev/stdout
             | cgat psl2psl
                  --queries-tsv-file=%(genome_query)s
                  --target-psl-file=%(genome_target)s
                  --method=sanitize
             | gzip
             >> %(outfile)s
             '''
            P.run()
else:
    raise ValueError(
        "configuration error: please specify either maf_dir or axt_dir")


def importRepeatsFromUCSC(infile, outfile, ucsc_database, repeattypes, genome):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses = "','".join(repeattypes.split(","))

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is
    # queried for tables that end in rmsk.
    dbhandle = PipelineUCSC.connectToUCSC(
        host=PARAMS["ucsc_host"],
        user=PARAMS["ucsc_user"],
        database=ucsc_database)

    cc = dbhandle.execute("SHOW TABLES LIKE '%%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    tmpfile = P.getTempFile(shared=True)

    total_repeats = 0
    for table in tables:
        E.info("%s: loading repeats from %s" % (ucsc_database, table))
        cc = dbhandle.execute(
            """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, '.',
            strand, '.',
            CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\";')
            FROM %(table)s
            WHERE repClass in ('%(repclasses)s') """ % locals())
        n = 0
        for data in cc.fetchall():
            n += 1
            tmpfile.write("\t".join(map(str, data)) + "\n")
        E.info("%s: %s=%i repeats downloaded" % (ucsc_database, table, n))
        total_repeats += n

    if total_repeats == 0:
        raise ValueErrror("did not find any repeats for %s" % ucsc_database)

    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''cat %(tmpfilename)s
    | %(pipeline_scriptsdir)s/gff_sort pos
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome)s
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()

    os.unlink(tmpfilename)


def importRepeatsFromEnsembl(infile, outfile,
                             ensembl_database,
                             repeattypes, genome):
    '''import repeats from an ENSEMBL database.
    '''
    statement = '''
    perl %(scriptsdir)s/ensembl_repeats2gff.pl
    -h %(ensembl_host)s
    -u %(ensembl_user)s
    -p %(ensembl_password)s
    -d %(ensembl_database)s
    --repeattypes %(repeattypes)s
    | %(pipeline_scriptsdir)s/gff_sort pos
    | cgat gff2gff
    --method=sanitize
    --sanitize-method=genome
    --skip-missing
    --genome-file=%(genome)s
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


@jobs_limit(1, "UCSC")
@files([(None, "%s_repeats.gff.gz" % x, x)
        for x in (PARAMS["query"], PARAMS["target"])])
def importRepeats(infile, outfile, track):
    '''import repeats from external sources.'''

    source = PARAMS["%s_source" % track]
    genome = os.path.join(PARAMS["genome_dir"], track)

    if source == "ensembl":
        importRepeatsFromEnsembl(infile, outfile,
                                 PARAMS["%s_database" % track],
                                 repeattypes=PARAMS["%s_repeattypes" % track],
                                 genome=genome)
    elif source == "ucsc":
        importRepeatsFromUCSC(infile, outfile,
                              PARAMS["%s_database" % track],
                              repeattypes=PARAMS["%s_repeattypes" % track],
                              genome=genome)


@transform(importRepeats,
           suffix("_repeats.gff.gz"),
           "_merged.gff.gz")
def mergeRepeats(infile, outfile):
    '''merge adjacent repeats.'''

    statement = '''gunzip
    < %(infile)s
    | cgat gff2gff
    --method=merge-features
    --min-distance=0
    --max-distance=10
    --min-features=0
    --max-features=0
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


@follows(buildGenomeAlignment)
@merge(mergeRepeats, "aligned_repeats.psl.gz")
def buildAlignedRepeats(infiles, outfile):
    '''build alignment between repeats.
    '''

    infile_target = PARAMS["target"] + "_merged.gff.gz"
    infile_query = PARAMS["query"] + "_merged.gff.gz"

    # using farm.py to send to cluster
    # granularity should be set automatically.
    granularity = 5000

    # need to escape pipe symbols within farm.py command
    # statement = r'''
    #     gunzip < %(interface_alignment_psl)s
    #     | %(cmd-farm)s --split-at-lines=%(granularity)i --log=%(outfile)s.log --is-binary
    #          "cgat psl2psl
    #             --method=test
    #     	--log=%(outfile)s.log
    #           | cgat psl2psl
    #     	--method=map
    #     	--filter-query=%(infile_query)s
    #     	--filter-target=%(infile_target)s
    #     	--log=%(outfile)s.log "
    #      | gzip
    #      > %(outfile)s'''
    # P.run()

    statement = '''
    gunzip < %(interface_alignment_psl)s
    | cgat psl2psl
    --method=test
    --log=%(outfile)s.log
    | cgat psl2psl
    --method=map
    --filter-query=%(infile_query)s
    --filter-target=%(infile_target)s
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s'''
    P.run()

########################################################
########################################################
########################################################


@files(buildAlignedRepeats, "aligned_repeats.rates.gz")
def buildRepeatsRates(infile, outfile):
    '''compute rates for individual aligned repeats.'''

    genome_query, genome_target = getGenomes()

    statement = '''gunzip < %(infile)s
    | sort -k10,10 -k14,14 -k9,9 -k12,12n
    | %(cmd-farm)s --split-at-lines=10000 --output-header --log=%(outfile)s.log
    "cgat psl2psl
    --log=%(outfile)s.log
    --method=add-sequence
    --queries-tsv-file=%(genome_query)s
    --target-psl-file=%(genome_target)s
    | cgat psl2table
    --method=query-counts
    --method=baseml
    --baseml-model=REV"
    | gzip > %(outfile)s
    '''
    P.run()


@transform((buildAlignedRepeats, buildGenomeAlignment),
           suffix(".psl.gz"),
           ".stats")
def computeAlignmentStats(infile, outfile):
    '''compute alignment coverage statistics'''

    statement = '''
    gunzip < %(infile)s
    | cgat psl2stats
    --log=%(outfile)s.log
    > %(outfile)s'''

    P.run()

########################################################
########################################################
########################################################


@transform(mergeRepeats, suffix(".gff.gz"), ".stats")
def computeRepeatsCounts(infile, outfile):
    '''count number and type of repeats.'''
    pass

# %_repeats_counts.stats: ucsc_%_repeats.table.gz
# 	$(PRELOG)
# @gunzip < $< | pe "s/#//" |\
# 	csv_cut genoName genoStart genoEnd repName repClass repFamily |\
# 	awk '/genoName/ {printf("%s\t%s\n", $$5, "length"); next;} {printf("%s\t%i\n", $$5, $$3-$$2); } ' |\
# 	t2t --group=1 --group-function=stats > $@
# 	$(EPILOG)

########################################################
########################################################
########################################################


@transform(mergeRepeats,
           suffix("_merged.gff.gz"),
           "_repeats_sizes.stats")
def buildRepeatDistribution(infile, outfile):
    '''count size and distance distribution of repeats.'''

    statement = '''gunzip
    < %(infile)s
    | cgat gff2histogram
    --output-filename-pattern="%(outfile)s.%%s"
    --method=all
    > %(outfile)s
    '''
    P.run()

########################################################
########################################################
########################################################


@files(buildRepeatsRates, PARAMS["interface_rates_query_gff"])
def exportRatesAsGFF(infile, outfile):
    '''export gff file with rate as score.'''

    statement = '''gunzip
    < %(infile)s
    | cgat csv_cut qName qStart qEnd distance converged
    | awk '!/qName/ && $5 {printf("%%s\\tancestral_repeat\\texon\\t%%s\\t%%s\\t%%s\\t+\\t.\\t.\\n", $1, $2, $3, $4);}'
    | gzip
    > %(outfile)s
    '''
    P.run()


@follows(importRepeats,
         mergeRepeats,
         buildAlignedRepeats,
         buildRepeatsRates,
         buildRepeatDistribution,
         computeAlignmentStats,
         computeRepeatsCounts,
         exportRatesAsGFF,
         )
def full():
    pass

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
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
