"""
pipeline_peakcalling.py - Window based genomic analysis
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python


Methods
=======

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----


Pipeline output
===============


Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *
import sys
import os
import re
import csv
import sqlite3
import glob
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelinePeakcallingScrum as PipelinePeakcalling
import CGAT.BamTools as Bamtools
import itertools
import pandas as pd
import numpy as np

#########################################################################
#########################################################################
#########################################################################
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
    update_interface=True))


# load IDR parameters into a dictionary to pass to the IDR step
# IDR requires multiple parameters from the PARAMS dictionary
idrPARAMS = dict()
idrpc = PARAMS['peakcalling_idrpeakcaller']
idrPARAMS['idrsuffix'] = PARAMS["%s_idrsuffix" % idrpc]
idrPARAMS['idrcol'] = PARAMS["%s_idrcol" % idrpc]
idrPARAMS['idrcolname'] = PARAMS['%s_idrcolname' % idrpc]
idrPARAMS['useoracle'] = PARAMS['IDR_useoracle']


# This function reads the design table and generates

# 1. A dictionary, inputD, linking each input file and each of the various
# IDR subfiles to the appropriate input, as specified in the design table

# 2. A pandas dataframe, df, containing the information from the
# design table
df, inputD = PipelinePeakcalling.readDesignTable("design.tsv",
                                                 PARAMS['IDR_poolinputs'])


# INPUTBAMS - list of control (input) bam files
# CHIPBAMS - list of experimental bam files on which to call peaks and perform
# IDR
INPUTBAMS = list(set(df['bamControl'].values))
CHIPBAMS = list(set(df['bamReads'].values))


# Check if reads are paired end
if Bamtools.isPaired(CHIPBAMS[0]) is True:
    PARAMS['paired_end'] = True
else:
    PARAMS['paired_end'] = False


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' %\
                (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###############################################################
# Preprocessing Steps


@follows(mkdir("filtered_bams.dir"))
@transform(INPUTBAMS, regex("(.*).bam"),
           [r"filtered_bams.dir/\1_filtered.bam",
            r"filtered_bams.dir/\1_counts.tsv"])
def filterInputBAMs(infile, outfiles):
    '''
    Applies various filters specified in the pipeline.ini to the bam file
    Currently implemented are filtering unmapped, unpaired and duplicate
    reads and filtering reads overlapping with blacklisted regions
    based on a list of bam files.
    '''
    filters = PARAMS['filters_bamfilters'].split(",")
    bedfiles = PARAMS['filters_bedfiles'].split(",")
    blthresh = PARAMS['filters_blacklistthresh']
    PipelinePeakcalling.filterBams(infile, outfiles, filters, bedfiles,
                                   float(blthresh),
                                   PARAMS['paired_end'],
                                   PARAMS['filters_strip'],
                                   PARAMS['filters_qual'],
                                   PARAMS['filters_keepint'])


@follows(mkdir("filtered_bams.dir"))
@transform(CHIPBAMS, regex("(.*).bam"), [r"filtered_bams.dir/\1_filtered.bam",
                                         r"filtered_bams.dir/\1_counts.tsv"])
def filterChipBAMs(infile, outfiles):
    '''
    Applies various filters specified in the pipeline.ini to the bam file
    Currently implemented are filtering unmapped, unpaired and duplicate
    reads and filtering reads overlapping with blacklisted regions
    based on a list of bam files.
    '''
    filters = PARAMS['filters_bamfilters'].split(",")
    bedfiles = PARAMS['filters_bedfiles'].split(",")
    blthresh = PARAMS['filters_blacklistthresh']
    PipelinePeakcalling.filterBams(infile, outfiles, filters, bedfiles,
                                   float(blthresh),
                                   PARAMS['paired_end'],
                                   PARAMS['filters_strip'],
                                   PARAMS['filters_qual'],
                                   PARAMS['filters_keepint'])


# These steps are required for IDR and are only run if IDR is requested
if int(PARAMS['IDR_run']) == 1:
    @follows(mkdir("pooled_bams.dir"))
    @split(filterChipBAMs,
           r"pooled_bams.dir/*_pooled_filtered.bam")
    def makePooledBams(infiles, outfiles):
        '''
        IDR requires one bam file for each replicate and a pooled bam
        file of all replicates for a particular condition and tissue.
        This function generates the pooled bam files.
        '''
        cond_tissues = set(df['Condition'] + "_" + df['Tissue'])

        # Take each combination of tissues and conditions from the design
        # tables
        for ct in cond_tissues:
            p = ct.split("_")
            cond = p[0]
            tissue = p[1].split(".")[0]

            # identify and read all bam files for this combination of
            # tissue and condition
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]
            innames = subdf['bamReads'].values
            innames = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in innames])

            out = "pooled_bams.dir/%s_pooled_filtered.bam" % ct

            # Generate a merged, sorted, indexed bam file combining
            # all bam files for this tissue and condition
            PipelinePeakcalling.mergeSortIndex(innames, out)

    @active_if(PARAMS['IDR_poolinputs'] != "all")
    @follows(mkdir('IDR_inputs.dir'))
    @split(filterInputBAMs, "IDR_inputs.dir/*_pooled_filtered.bam")
    def makePooledInputs(infiles, outfiles):
        '''
        As pooled BAM files are used in the IDR, pooled input files also
        need to be generated - combined bam files of all the input bam
        files for this tissue.
        If you have chosen the "all" option for IDR_poolinputs in the
        pipeline.ini, this step is skipped, as all inputs are pooled for
        all IDR analyses.
        '''
        cond_tissues = set(df['Condition'] + "_" + df['Tissue'])

        # Take each combination of tissues and conditions from the design
        # tables
        for ct in cond_tissues:
            p = ct.split("_")
            cond = p[0]
            tissue = p[1].split(".")[0]
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]

            # find the inputs linked to any bam files for this combination of
            # tissues and conditions
            inputs = subdf['bamControl'].values
            inputs = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in inputs])
            out = "IDR_inputs.dir/%s_pooled_filtered.bam" % ct

            # generate a sorted, index, merged bam file for all of these
            # inputs
            PipelinePeakcalling.mergeSortIndex(inputs, out)

else:
    @transform(filterChipBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'filtered_bams.dir/\1.bam')
    def makePooledBams(infile, outfile):
        '''
        Dummy task if IDR not requested.
        '''
        pass

    @transform(filterInputBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'filtered_bams.dir/\1.bam')
    def makePooledInputs(infile, outfile):
        pass


if int(PARAMS['IDR_run']) == 1:
    @follows(mkdir("peakcalling_bams.dir"))
    @subdivide((filterChipBAMs, makePooledBams),
               regex("(.*)_bams.dir/(.*).bam"),
               [r"peakcalling_bams.dir/\2_pseudo_1.bam",
                r"peakcalling_bams.dir/\2_pseudo_2.bam",
                r"peakcalling_bams.dir/\2.bam"])
    def makePseudoBams(infiles, outfiles):
        '''
        Generates pseudo bam files each containing approximately 50% of reads
        from the original bam file for IDR self consistency analysis.
        Also generates a link to the original BAM file in the
        peakcalling_bams.dir directory.

        '''
        # makePooledBams generates a single output whereas filterChipBAMS
        # generates a bam file and a table - a list of outputs
        if isinstance(infiles, list):
            infile = infiles[0]
        else:
            infile = infiles

        pseudos = outfiles[0:2]
        orig = outfiles[2]

        PipelinePeakcalling.makeBamLink(infile, orig)

        PipelinePeakcalling.makePseudoBams(infile, pseudos,
                                           PARAMS['paired_end'],
                                           PARAMS['IDR_randomseed'],
                                           PARAMS['filters_bamfilters'].split(
                                               ","),
                                           submit=True)
else:
    @follows(mkdir('peakcalling_bams.dir'))
    @transform(filterChipBAMs, regex("filtered_bams.dir/(.*)_filtered.bam"),
               r'peakcalling_bams.dir/\1.bam')
    def makePseudoBams(infile, outfile):
        '''
        Link to original BAMs without generating pseudo bams
        if IDR not requested.
        '''
        PipelinePeakcalling.makeBamLink(infile[0], outfile)


# These three functions gather and parse the input (control) bam files into the
# IDR_inputs.dir directory prior to IDR analysis.
# The method used to do this depends on the IDR_poolinputs parameter

if PARAMS['IDR_poolinputs'] == "none":
    @follows(mkdir('IDR_inputs.dir'))
    @transform(filterInputBAMs, regex("filtered_bams.dir/(.*).bam"),
               r'IDR_inputs.dir/\1.bam')
    def makeIDRInputBams(infile, outfile):
        '''
        When pooled inputs are not requested, the appropriate inputs are
        generated above in the filterInputBAMS step - this function links to
        these in the IDR_inputs.dir directory.
        '''
        infile = infile[0]
        PipelinePeakcalling.makeBamLink(infile, outfile)


elif PARAMS['IDR_poolinputs'] == "all":
    @follows(mkdir('IDR_inputs.dir'))
    @merge(filterInputBAMs, "IDR_inputs.dir/pooled_all.bam")
    def makeIDRInputBams(infiles, outfile):
        '''
        When all inputs are to be pooled and used as a control against all
        samples, a single merged bam is generated from the output of
        the filterInputBAMs step above in the IDR_inputs.dir directory.
        '''
        infiles = [i[0] for i in infiles]
        PipelinePeakcalling.mergeSortIndex(infiles, outfile)


elif PARAMS['IDR_poolinputs'] == "condition" and PARAMS['IDR_run'] != 1:
    @follows(mkdir('IDR_inputs.dir'))
    @split(filterInputBAMs, r'IDR_inputs.dir/*.bam')
    def makeIDRInputBams(infiles, outfiles):
        '''
        When IDR is going to be performed, inputs which are pooled by tissue
        and condition are automatically generated as these are always required.

        This function pools tissues and conditions when IDR is switched
        off if inputs pooled by condition are requested.

        The appropriate outputs from filterInputBAMs are identified and
        merged into a single BAM stored in the IDR_inputs.dir directory.
        '''
        outs = set(inputD.values())
        for out in outs:
            p = out.split("_")
            cond = p[1]
            tissue = p[2].split(".")[0]

            # collect the appropriate bam files from their current location
            subdf = df[((df['Condition'] == cond) & (df['Tissue'] == tissue))]
            innames = subdf['bamControl'].values
            innames = set(
                ["filtered_bams.dir/%s" % s.replace(".bam", "_filtered.bam")
                 for s in innames])
            out = "IDR_inputs.dir/%s" % out

            PipelinePeakcalling.mergeSortIndex(innames, out)


elif PARAMS['IDR_poolinputs'] == "condition" and PARAMS['IDR_run'] == 1:
    @follows(mkdir('IDR_inputs.dir'))
    @follows(mkdir('IDR_inputs.dir'))
    @transform(makePooledInputs, regex("IDR_inputs.dir/(.*).bam"),
               r'IDR_inputs.dir/\1.bam')
    def makeIDRInputBams(infiles, outfiles):
        '''
        If IDR is going to be run, pooled inputs are generated above so
        they don't need to be generated again if requested.
        '''
        pass


@follows(makeIDRInputBams)
@follows(filterInputBAMs)
@follows(makePooledBams)
@follows(makePooledInputs)
@follows(makePseudoBams)
@originate("peakcalling_bams_and_inputs.tsv")
def makeBamInputTable(outfile):
    '''
    Generates a tab delimited file - peakcalling_bams_and_inputs.tsv
    which links each filtered bam file in the peakcalling_bams.dir
    directory to the appropriate input in the IDR_inputs.dir
    directory.
    Uses the dictionary inputD generated as a global variable based
    on the user-specified design table plus pooled input files generated
    above.
    '''
    ks = inputD.keys()
    out = IOTools.openFile(outfile, "w")
    out.write('ChipBam\tInputBam\n')
    bamfiles = os.listdir("peakcalling_bams.dir")

    for k in ks:
        inputstem = inputD[k]
        chipstem = k
        chipstem = P.snip(chipstem)
        inputstem = P.snip(inputstem)
        inputfile = "IDR_inputs.dir/%s_filtered.bam" % inputstem

        for b in bamfiles:
            if b.startswith(chipstem) and b.endswith('bam'):
                out.write("peakcalling_bams.dir/%s\t%s\n" % (b, inputfile))
    out.close()


@transform(makePseudoBams, suffix(".bam"), "_insertsize.tsv")
def estimateInsertSize(infile, outfile):
    '''
    Predicts insert size using MACS2 for single end data and using Bamtools
    for paired end data.
    Output is stored in insert_size.tsv
    '''
    PipelinePeakcalling.estimateInsertSize(infile, outfile,
                                           PARAMS['paired_end'],
                                           PARAMS['insert_alignments'],
                                           PARAMS['insert_macs2opts'])


@merge(estimateInsertSize, "insert_sizes.tsv")
def mergeInsertSizes(infiles, outfile):
    '''
    Combines insert size outputs into one file
    '''
    out = IOTools.openFile(outfile, "w")
    out.write("filename\tmode\tfragmentsize_mean\tfragmentsize_std\ttagsize\n")
    for infile in infiles:
        res = IOTools.openFile(infile).readlines()
        out.write("%s\t%s\n" % (infile, res[-1].strip()))
    out.close()


@follows(makeBamInputTable)
@follows(mergeInsertSizes)
@transform(makePseudoBams, regex("(.*)_bams\.dir\/(.*)\.bam"),
           r"\1_bams.dir/\2.bam")
def preprocessing(infile, outfile):
    '''
    Dummy task to ensure all preprocessing has run and
    bam files are passed individually to the next stage.
    '''
    pass

# ###############################################################
#  Peakcalling Steps


@follows(mkdir('macs2.dir'))
@transform(preprocessing,
           regex("peakcalling_bams.dir/(.*).bam"),
           add_inputs(makeBamInputTable),
           r"macs2.dir/\1.macs2")
def callMacs2peaks(infiles, outfile):
    D = PipelinePeakcalling.readTable(infiles[1])
    bam = infiles[0]
    inputf = D[bam]
    insertsizef = "%s_insertsize.tsv" % (P.snip(bam))

    peakcaller = PipelinePeakcalling.Macs2Peakcaller(
        threads=1,
        paired_end=True,
        output_all=True,
        tool_options=PARAMS['macs2_options'],
        tagsize=None)

    statement = peakcaller.build(bam, outfile,
                                 PARAMS['annotations_interface_contigs'],
                                 inputf, insertsizef, PARAMS['IDR_run'],
                                 PARAMS['macs2_idrkeeppeaks'],
                                 PARAMS['macs2_idrsuffix'],
                                 PARAMS['macs2_idrcol'])
    P.run()


PEAKCALLERS = []
IDRPEAKCALLERS = []
mapToPeakCallers = {'macs2': (callMacs2peaks,)}

for x in P.asList(PARAMS['peakcalling_peakcallers']):
    PEAKCALLERS.extend(mapToPeakCallers[x])


@follows(*PEAKCALLERS)
def peakcalling():
    '''
    dummy task to define upstream peakcalling tasks
    '''

################################################################
# IDR Steps

@follows(peakcalling)
@follows(mkdir("peaks_for_IDR.dir"))
@transform(mapToPeakCallers[PARAMS['peakcalling_idrpeakcaller']],
           regex("(.*)/(.*)"),
           r"peaks_for_IDR.dir/\2.IDRpeaks")
def getIDRInputs(infile, outfile):
    IDRpeaks = "%s_IDRpeaks" % infile
    shutil.copy(IDRpeaks, outfile)


@merge(getIDRInputs, "IDR_pairs.tsv")
def makeIDRPairs(infiles, outfile):
    useoracle = PARAMS['IDR_useoracle']
    PipelinePeakcalling.makePairsForIDR(infiles, outfile,
                                        PARAMS['IDR_useoracle'],
                                        df, submit=True)

@follows(mkdir("IDR.dir"))
@split(makeIDRPairs, "IDR.dir/*.dummy")
def splitForIDR(infile, outfiles):
    pairs = pd.read_csv(infile, sep="\t", header=None)
    for p in pairs.index.values:
        p = pairs.ix[p]
        p1 = P.snip(p[0].split("/")[-1])
        p2 = P.snip(p[1].split("/")[-1])

        pairstring = "%s_v_%s" % (p1, p2)

        out = IOTools.openFile("IDR.dir/%s.dummy" % pairstring, "w")
        out.write("%s\n" % "\n".join(p))
        out.close()


@transform(splitForIDR, suffix(".dummy"), ".tsv")
def runIDR(infile, outfile):
    lines = [line.strip() for line in IOTools.openFile(infile).readlines()]
    infile1, infile2, setting, oraclefile, condition, tissue = lines
    options = PARAMS['IDR_options']

    if setting == 'self_consistency':
        idrthresh = PARAMS['IDR_softthresh_selfconsistency']
        options += " %s" % PARAMS['IDR_options_selfconsistency']
    elif setting == "pooled_consistency":
        idrthresh = PARAMS['IDR_softthresh_pooledconsistency']
        options += " %s" % PARAMS['IDR_options_pooledconsistency']

    elif setting == "replicate_consistency":
        idrthresh = PARAMS['IDR_softthresh_replicateconsistency']
        options += " %s" % PARAMS['IDR_options_replicateconsistency']

    statement = PipelinePeakcalling.buildIDRStatement(
        infile1, infile2,
        outfile,
        PARAMS['IDR_sourcecommand'],
        PARAMS['IDR_unsourcecommand'],
        idrthresh,
        idrPARAMS, options, oraclefile)

    P.run()


@transform(runIDR, suffix(".tsv"), ["_filtered.tsv",
                                    "_table.tsv"])
def filterIDR(infile, outfiles):
    '''
    Take the IDR output, which is in ENCODE narrowPeaks format if the input
    is narrowPeaks, gtf or bed and ENCODE broadPeaks format if the input is
    broadPeaks.
    Input is filtered based on whether it passes the soft IDR thresholds
    provided in the pipeline.ini.  Peaks which pass this threshold
    with have a score in the "globalIDR" column which is greater
    than -log(soft_threshold) where soft_threshold is the soft threshold
    provided in the pipeline.ini.
    Column headings are added and output is sorted by signalValue.
    '''
    IDRdata = pd.read_csv(infile, sep="\t", header=None)
    if idrPARAMS['idrsuffix'] == "broadPeaks":
        IDRdata.columns = ["chrom", "chromStart", "chromEnd", "name", "score",
                           "strand", "signalValue", "p-value", "q-value",
                           "localIDR", "globalIDR",
                           "rep1_chromStart", "rep2_chromEnd",
                           "rep1_signalValue", "rep2_chromStart",
                           "rep2_chromEnd", "rep2_signalValue"]
    else:
        IDRdata.columns = ["chrom", "chromStart", "chromEnd", "name", "score",
                           "strand", "signalValue", "p-value", "q-value",
                           "summit", "localIDR", "globalIDR",
                           "rep1_chromStart", "rep2_chromEnd",
                           "rep1_signalValue", "rep1_summit",
                           "rep2_chromStart", "rep2_chromEnd",
                           "rep2_signalValue", "rep2_summit"]

    IDRdataP = IDRdata[IDRdata['score'] == 1000]
    IDRdataF = IDRdata[IDRdata['score'] != 1000]

    IDRdataP = IDRdataP.sort_values('signalValue', ascending=False)
    IDRdataF = IDRdataF.sort_values('signalValue', ascending=False)

    H = ['Total_Peaks', 'Peaks_Passing_IDR', 'Peaks_Failing_IDR',
         'Percentage_Peaks_Passing_IDR']
    T = ((len(IDRdata), len(IDRdataP), len(IDRdataF),
          round(float(len(IDRdataP)) / float(len(IDRdata)), 4) * 100))

    out = IOTools.openFile(outfiles[1], "w")
    out.write("%s\n" % "\t".join(H))
    out.write("%s\n" % "\t".join([str(t) for t in T]))

    IDRdataP.to_csv(outfiles[0], sep="\t")


@merge((filterIDR, makeIDRPairs), "IDR_results.tsv")
def summariseIDR(infiles, outfile):
    tables = [i[1] for i in infiles[:-1]]
    pairs = infiles[-1]

    alltab = pd.DataFrame()
    for table in tables:
        tab = pd.read_csv(table, sep="\t")
        tab['Output_Filename'] = table.split("/")[-1]
        alltab = alltab.append(tab)
    pairtab = pd.read_csv(pairs, sep="\t", names=['Input_File_1',
                                                  'Input_File_2',
                                                  'Replicate_Type',
                                                  'Oracle_Peak_File',
                                                  'Condition',
                                                  'Tissue'])

    # Generate a temporary column to merge the two tables
    pairstrings = []
    for p in pairtab.index.values:
        p = pairtab.ix[p]
        p1 = P.snip(p[0].split("/")[-1])
        p2 = P.snip(p[1].split("/")[-1])
        pairstring = "%s_v_%s_table.tsv" % (p1, p2)
        pairstrings.append(pairstring)

    pairtab['pairstring'] = pairstrings
    alltab = alltab.merge(pairtab, left_on='Output_Filename',
                          right_on='pairstring')
    alltab = alltab.drop('pairstring', 1)

    thresholds = []
    for item in alltab['Replicate_Type']:
        if item == "pooled_consistency":
            thresholds.append(PARAMS['IDR_softthresh_pooledconsistency'])
        elif item == "self_consistency":
            thresholds.append(PARAMS['IDR_softthresh_selfconsistency'])
        elif item == "replicate_consistency":
            thresholds.append(PARAMS['IDR_softthresh_replicateconsistency'])
    print thresholds
    alltab['IDR_Thresholds'] = thresholds
    alltab['IDR_Thresholds'] = alltab['IDR_Thresholds'].astype('float')
    alltab['IDR_Transformed_Thresholds'] = -(np.log10(thresholds))

    alltab = alltab[['Output_Filename', 'Replicate_Type', 'Total_Peaks',
                     'Peaks_Passing_IDR', 'Percentage_Peaks_Passing_IDR',
                     'Peaks_Failing_IDR', 'IDR_Thresholds',
                     'IDR_Transformed_Thresholds',
                     'Input_File_1', 'Input_File_2', 'Condition', 'Tissue',
                     'Oracle_Peak_File']]
    alltab.to_csv(outfile, sep="\t", index=False)

@transform(summariseIDR, suffix("results.tsv"), "analysis.tsv")
def analyseIDR(infile, outfile):
    table = pd.read_csv(infile, sep="\t")
    

################################################################
# QC Steps


################################################################
# Fragment GC% distribution
################################################################

"""
@follows(mkdir("QC.dir"))
@transform(BAMS, regex("(.*).bam"), r"QC.dir/\1.tsv")
def fragLenDist(infile, outfile):

    if PARAMS["paired_end"] == 1:
        function = "--merge-pairs"
    else:
        function = "--fragment"

    genome = os.path.join(PARAMS["general_genome_dir"],
                          PARAMS["general_genome"])
    genome = genome + ".fasta"

    statement = '''
    samtools view -s 0.2 -ub %(infile)s |
    python %(scriptsdir)s/bam2bed.py  %(function)s |
    python %(scriptsdir)s/bed2table.py --counter=composition-na -g %(genome)s\
    > %(outfile)s
    '''
    P.run()


@merge(BAMS, regex("(.*).bam"), r"QC.dir/genomic_coverage.tsv")
def buildReferenceNAComposition(infiles, outfile):

    infile = infiles[0]
    contig_sizes = os.path.join(PARAMS["annotations_dir"],
                                PARAMS["annotations_interface_contigs"])
    gaps_bed = os.path.join(PARAMS["annotations_dir"].
                            PARAMS["annotations_interface_gaps_bed"])

    statement = '''bedtools shuffle
    -i %(infile)s
    -g %(contig_sizes)s
    -excl %(gaps_bed)s
    -chromFirst
    | python %(scriptsdir)s/bed2table.py
    --counter=composition-na
    -g %(genome)s > %(outfile)s
    '''

    P.run()
################################################################
"""


def full():
    pass


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


@follows(mkdir("%s/bamfiles" % PARAMS["web_dir"]),
         mkdir("%s/medips" % PARAMS["web_dir"]),
         )
def publish():
    '''publish files.'''

    # directory : files

    # publish web pages
    P.publish_report(export_files=export_files)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
