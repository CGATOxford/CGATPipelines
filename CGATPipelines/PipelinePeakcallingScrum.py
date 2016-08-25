'''
======================================================================

Requirements:


Reference
---------

'''

import os
import re
import collections
import itertools
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools
import pandas as pd
import pysam
import numpy as np
import shutil
from CGATPipelines.Pipeline import cluster_runnable
import rpy2
from rpy2.robjects import r as R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
pandas2ri.activate()


##############################################
# Preprocessing Functions


def readDesignTable(infile, poolinputs):
    '''
    This function reads the design table and generates

    1. A dictionary, inputD, linking each input file and each of the various
    IDR subfiles to the appropriate input, as specified in the design table

    2. A pandas dataframe, df, containing the information from the
    design table

    This dictionary includes the filenames that the IDR output files will have
    once these steps of the pipeline are run
    If IDR is not requsted, these dictionary values will still exist but will
    not be used.

    poolinputs has three possible options:

    none
    Use the input files per replicate as specified in the design files
    Input files will still also be pooled if IDR is specified as IDR requires
    BAM files representing pooled replicates as well as BAM files for each
    replicate

    all
    Pool all input files and use this single pooled input BAM file as the input
    for any peakcalling.  Used when input is low depth or when only a single
    input or replicates of a single input are available.

    condition
    pool the input files within each combination of conditions and tissues
    '''

    df = pd.read_csv("design.tsv", sep="\t")
    pairs = zip(df['bamReads'], df['bamControl'])
    CHIPBAMS = list(set(df['bamReads'].values))

    if poolinputs == "none":
        inputD = dict(pairs)
        conditions = df['Condition'].values
        tissues = df['Tissue'].values
        i = 0
        for C in CHIPBAMS:
            cond = conditions[i]
            tissue = tissues[i]
            inputD["%s_%s.bam" % (cond, tissue)] = "%s_%s_pooled.bam" % (
                cond, tissue)
            i += 1

    elif poolinputs == "all":
        inputD = dict()
        for C in CHIPBAMS:
            inputD[C] = "pooled_all.bam"

    elif poolinputs == "condition":
        inputD = dict()
        conditions = df['Condition'].values
        tissues = df['Tissue'].values
        i = 0
        for C in CHIPBAMS:
            cond = conditions[i]
            tissue = tissues[i]
            inputD[C] = "%s_%s_pooled.bam" % (cond, tissue)
            inputD["%s_%s.bam" % (cond, tissue)] = "%s_%s_pooled.bam" % (
                cond, tissue)
            i += 1
    df['Condition'] = df['Condition'].astype('str')
    df['Tissue'] = df['Tissue'].astype('str')
    return df, inputD


def trackFilters(filtername, bamfile, tabout):
    '''
    Keeps track of which filters are being used on the bam files and generates
    a statement segment which records this in a log file
    '''
    return """echo %(filtername)s >> %(tabout)s;
              samtools view -c %(bamfile)s.bam >> %(tabout)s; """ % locals()


def appendSamtoolsFilters(statement, inT, tabout, filters, qual, pe):
    '''
    Apply filters using samtools
    -F - remove these flags
    -f - keep only these flags
    -f1, -f2 remove reads which are unpaired or not properly paired
    -F4 removes unmapped reads
    '''
    i = 0
    j = 0
    string = ""
    for filt in filters:
        if filt == "unmapped":
            string = "-F 4"
        elif filt == "secondary":
            string = "-F 0x100"
        elif filt == "unpaired":
            if pe is True:
                string = "-f 1 -f 2"
            else:
                E.warn("""Bam file is not paired-end so unpaired reads have
                not been filtered""")
                continue
        elif filt == "lowqual":
            string = '-q %(qual)i' % locals()
        # append any new samtools filters to this list
        if filt in ["unmapped", "unpaired", "lowqual", "secondary"]:
            if i == 0:
                outT = P.getTempFilename("./filtered_bams.dir")
            else:
                inT = outT
                outT = P.getTempFilename("./filtered_bams.dir")

            statement += """samtools view -b %(string)s %(inT)s.bam
            > %(outT)s.bam; rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()
            statement += trackFilters(filt, outT, tabout)
            i += 1

    statement = statement.replace("\n", "")
    return statement, outT


def appendPicardFilters(statement, inT, tabout, filters, pe, outfile):
    '''
    Apply filters using Picard.
    MarkDuplicates removes duplicate reads.
    '''
    outT = inT
    if 'duplicates' in filters:
        log = outfile.replace(".bam", "_duplicates.log")
        outT = P.getTempFilename("./filtered_bams.dir")
        statement += """
        MarkDuplicates \
        INPUT=%(inT)s.bam \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        OUTPUT=%(outT)s.bam \
        METRICS_FILE=/dev/null \
        VALIDATION_STRINGENCY=SILENT \
        2> %(log)s; rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()

        statement += trackFilters("duplicates", outT, tabout)
        statement = statement.replace("\n", "")
    return statement, outT


def appendBlacklistFilter(statement, inT, tabout, bedfiles, blthresh, pe):
    '''
    Filters blacklisted regions based on a list of bed files
    By default the read is removed if it has any overlap with the blacklisted
    region on either end.
    If the read is paired end both halves of the pair are removed if one
    overlaps with a blacklisted region
    '''
    outT = P.getTempFilename("./filtered_bams.dir")
    if pe is True:
        statement += """samtools sort -n %(inT)s.bam -o %(outT)s.bam;
                        rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()
        for bedfile in bedfiles:
            inT = outT
            outT = P.getTempFilename("./filtered_bams.dir")
            statement += """pairToBed -abam %(inT)s.bam
                                      -b %(bedfile)s
                                      -f %(blthresh)f10 -type neither
                            > %(outT)s.bam;
                            rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()
            statement += trackFilters(bedfile, outT, tabout)
            inT = outT
        outT = P.getTempFilename("./filtered_bams.dir")
        statement += """samtools sort %(inT)s.bam -o %(outT)s.bam;
                        rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()

    else:
        for bedfile in bedfiles:
            outT = P.getTempFilename("./filtered_bams.dir")
            statement += """bedtools intersect -abam %(inT)s.bam
                                      -b %(bedfile)s
                                      -f %(blthresh)f10 -v
                            > %(outT)s.bam;
                            rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()
            statement += trackFilters(bedfile, outT, tabout)
            inT = outT

    return statement, outT


def filterBams(infile, outfiles, filters, bedfiles, blthresh, pe, strip, qual,
               keep_intermediates=False):
    '''
    Applies various filters to bam files.
    Builds a statement to process the bam file - the file is sorted
    then the samtools and picard filters specified in the pipeline.ini
    are applied.  Reads overlapping with regions in the specified bed
    files are filtered out.  The filtered bam file is then indexed.  Read
    counts after each filtering step are stored in filtered_bams/*_counts.tsv.
    '''
    bamout, tabout = outfiles
    o = open(tabout, "w")
    o.close()
    inT = infile

    outT = P.getTempFilename("./filtered_bams.dir")

    statement = """samtools sort %(inT)s -o %(outT)s.bam; """ % locals()
    inT = outT

    statement += trackFilters("none", inT, tabout)

    statement, inT = appendPicardFilters(statement, inT, tabout, filters, pe,
                                         bamout)
    statement, inT = appendSamtoolsFilters(statement, inT, tabout, filters,
                                           qual, pe)
    statement, inT = appendBlacklistFilter(statement, inT, tabout, bedfiles,
                                           blthresh, pe)

    # I added isStripped to BamTools - this is commented out until it
    # is merged (Katy)
    # if int(strip) == 1 and BamTools.isStripped(inT) is False:
    #     # strip sequence if requested
    #     outT = P.getTempFilename(".")
    #     statement += """python %%(scriptsdir)s/bam2bam.py
    #                     -I %(inT)s.bam
    #                     --strip-method=all
    #                     --method=strip-sequence
    #                     --log=%(bamout)s.log -S %(outT)s.bam;
    #                     rm -f %(inT)s; """ % locals()

    #     inT = outT
    statement += """mv %(inT)s.bam %(bamout)s;
                    rm -f %(inT)s;
                    samtools index %(bamout)s""" % locals()

    statement = statement.replace("\n", "")

    if int(keep_intermediates) == 1:
        statement = re.sub("rm -f \S+.bam;", "", statement)
    P.run()

    # reformats the read counts into a table
    inf = [line.strip() for line in open(tabout).readlines()]
    i = 0
    ns = []
    nams = []
    for line in inf:
        line = line.strip()
        if len(line) != 0:
            if i % 2 == 0:
                ns.append(line)
            else:
                nams.append(line)
            i += 1
    o = open(tabout, "w")
    o.write("%s\n%s" % ("\t".join(ns), "\t".join(nams)))
    o.close()

    # check the filtering is done correctly - write a log file
    # if unpaired is specified in bamfilters in the pipeline.ini
    # remove reads whose mate has been filtered out elsewhere

    T = P.getTempFilename(".")
    checkBams(bamout, filters, qual, pe, T)
    if int(keep_intermediates) == 1:
        shutil.copy(bamout, bamout.replace(".bam", "_beforepaircheck.bam"))
    shutil.move("%s.bam" % T, bamout)
    shutil.move("%s.filteringlog" % (T),
                bamout.replace(".bam", ".filteringlog"))
    if os.path.exists(T):
        os.remove(T)
    sortIndex(bamout)


@cluster_runnable
def checkBams(infile, filters, qlim, pe, outfile):

    samfile = pysam.AlignmentFile(infile, 'rb')

    outbam = pysam.AlignmentFile("%s.bam" % outfile, "wb",
                                 template=samfile)
    sep = '.'
    logfile = sep.join([outfile, 'filteringlog'])

    counter = collections.Counter()
    fragment_length = collections.Counter()
    d = collections.defaultdict(list)

    if "lowqual" not in filters:
        qlim = 0

    counter['secondary'] = 0
    counter['proper_pairs'] = 0
    counter['improper_pairs'] = 0
    counter['high_quality'] = 0
    counter['low_quality'] = 0
    counter['unmapped'] = 0
    counter['mapped'] = 0
    counter['multiple_or_1_read_in_pair'] = 0

    for read in samfile.fetch():
        d[read.query_name].append(read)

        counter['total_reads'] += 1
        if read.is_secondary:
            counter['secondary'] += 1

        if read.is_proper_pair:
            counter['proper_pairs'] += 1
        else:
            counter['improper_pairs'] += 1

        if read.mapping_quality >= qlim:
            counter['high_quality'] += 1
        else:
            counter['low_quality'] += 1

        if read.is_unmapped:
            counter['unmapped'] += 1
        else:
            counter['mapped'] += 1

    if "unpaired" not in filters or pe is False:
        shutil.copy(infile, outfile)

    if "unpaired" in filters:
        for items, values in d.iteritems():
            if len(values) == 2:
                if values[0].is_read1 and values[1].is_read2:
                    outbam.write(values[0])
                    outbam.write(values[1])

                elif values[0].is_read2 and values[1].is_read1:
                    outbam.write(values[1])
                    outbam.write(values[0])

                else:
                    counter['paired11_paired22'] += 1
                l = abs(values[0].template_length)
                fragment_length[l] += 1

            else:
                counter['multiple_or_1_read_in_pair'] += 1
                if "secondary" not in filters:
                    for v in values:
                        outbam.write(v)

    out = IOTools.openFile(infile.replace(".bam", ".fraglengths"), "w")
    for key in fragment_length:
        out.write("%s\t%d\n" % (key, fragment_length[key]))
    out.close()

    filteringReport(counter, logfile)

    outbam.close()
    samfile.close()


def filteringReport(counter, logfile):
    logfile = open(logfile, "w")
    for c in counter:
        logfile.write("%s\t%s\n" % (c, counter[c]))
    logfile.close()


def estimateInsertSize(infile, outfile, pe, nalignments, m2opts):
    '''predict fragment size.

    For single end data, use MACS2, for paired-end data
    use the bamfile.

    In some BAM files the read length (rlen) attribute
    is not set, which causes MACS2 to predict a 0.

    Thus it is computed here from the CIGAR string if rlen is 0.
    '''

    print infile
    tagsize = BamTools.estimateTagSize(infile, multiple="mean")

    if pe is True:
        mode = "PE"
        mean, std, n = BamTools.estimateInsertSizeDistribution(
            infile, int(nalignments))
    else:
        mode = "SE"
        statement = '''macs2 predictd
        --format BAM
        --ifile %(infile)s
        --outdir %(outfile)s.dir
        --verbose 2
        %(insert_macs2opts)s
        >& %(outfile)s.tsv
        '''
        P.run()

        with IOTools.openFile(outfile + ".tsv") as inf:
            lines = inf.readlines()
            line = [x for x in lines
                    if "# predicted fragment length is" in x]
            if len(line) == 0:
                raise ValueError(
                    'could not find predicted fragment length')
            mean = re.search(
                "# predicted fragment length is (\d+)",
                line[0]).groups()[0]
            std = 'na'

    outf = IOTools.openFile(outfile, "w")
    outf.write("mode\tfragmentsize_mean\tfragmentsize_std\ttagsize\n")
    outf.write("\t".join(
        map(str, (mode, mean, std, tagsize))) + "\n")
    outf.close()


@cluster_runnable
def makePseudoBams(infile, outfiles, pe, randomseed, filters):
    '''
    Generates pseudo bam files for IDR analysis
    Each read in the input bam is assigned to a pseudo bam at random.
    If reads are paired end the pair will go to the same bam.
    '''

    # read bam file
    bamfile = pysam.AlignmentFile(infile, "rb")
    T = P.getTempFilename(".")
    # sort
    # pysam.sort("-n", infile, T, catch_stdout=False)
    os.system("""samtools sort -n %(infile)s -o %(T)s.bam""" % locals())

    sorted_bamfile = pysam.AlignmentFile("%s.bam" % T, "rb")

    # for single end, count the reads, for paired end, halve number of reads
    # then generate a random list of 0s and 1s of this length
    # 0 = go to pseudo bam 0, 1 = go to pseudo bam 1
    bamlength = bamfile.count()
    if pe is True:
        countreads = bamlength / 2
    else:
        countreads = bamlength
    randomgen = np.random.RandomState()
    randomgen.seed(randomseed)
    intlist = randomgen.random_integers(0, 1, countreads)
    outs = [pysam.AlignmentFile(outfiles[0], "wb", template=bamfile),
            pysam.AlignmentFile(outfiles[1], "wb", template=bamfile)]
    j = 0
    i = 0
    k = 0

    # send reads to replicates according to the integers in intlist
    dest = outfiles[i]
    if pe is True:
        for read in sorted_bamfile:
            # if j is even
            if j % 2 == 0:
                # take item i from intlist
                destint = intlist[i]
                # sends to output bam file 0 or 1
                dest = outs[destint]
                i += 1
            dest.write(read)
            j += 1
    else:
        for read in sorted_bamfile:
            destint = intlist[k]
            dest = outs[destint]
            dest.write(read)
            k += 1

    outs[0].close()
    outs[1].close()
    os.remove(T)
    os.remove("%s.bam" % T)
    lens = []

    # check that there are twice as many reads as read names for a paired
    # end bam file and check the lengths of the files
    for outf in outfiles:
        T = P.getTempFilename(".")
        os.system("""samtools sort %(outf)s -o %(T)s.bam;
        samtools index %(T)s.bam""" % locals())
  #      pysam.sort(outf, T, catch_stdout=False)
  #      pysam.index("%s.bam" % T, catch_stdout=False)

        bamfile = pysam.AlignmentFile("%s.bam" % T, "rb")
        all = []
        for nam in bamfile:
            all.append(nam.qname)

        if pe is True and "unpaired" in filters and "secondary" in filters:
            allreads = len(all)
            uniquereads = len(set(all))
            expectedreads = len(all) / 2
            assert (
                (len(set(all)) <= (len(all) / 2) + 2) &
                (len(set(all)) >= (len(all) / 2) - 2)), """
                Error splitting bam file %(outf)s\
                %(allreads)i reads in bam file and\
                %(uniquereads)s unique read names -
                expecting %(expectedreads)i unique read names\
                """ % locals()

        lens.append(bamfile.count())
        os.remove(T)
        shutil.move("%s.bam" % T, outf)
        shutil.move("%s.bam.bai" % T, "%s.bai" % outf)

    E.info("Bamfile 1 length %i, Bamfile 2 length %i" % (lens[0], lens[1]))

#############################################
# Peakcalling Functions


def getMacsPeakShiftEstimate(infile):
    ''' parse the peak shift estimate file from macs,
    return the fragment size '''

    with IOTools.openFile(infile, "r") as inf:

        header = inf.next().strip().split("\t")
        values = inf.next().strip().split("\t")

        fragment_size_mean_ix = header.index("fragmentsize_mean")

        fragment_size = int(float(values[fragment_size_mean_ix]))

        return fragment_size


def mergeSortIndex(bamfiles, out):
    '''
    Merge bamfiles into a single sorted, indexed outfile.
    '''
    infiles = " ".join(bamfiles)
    T1 = P.getTempFilename(".")
    T2 = P.getTempFilename(".")
    statement = """samtools merge %(T1)s.bam %(infiles)s;
    samtools sort %(T1)s.bam -o %(T2)s.bam;
    samtools index %(T2)s.bam;
    mv %(T2)s.bam %(out)s;
    mv %(T2)s.bam.bai %(out)s.bai""" % locals()
    P.run()
    os.remove("%s.bam" % T1)
    os.remove(T1)


def sortIndex(bamfile):
    '''
    Merge bamfiles into a single sorted, indexed outfile.
    '''
    T1 = P.getTempFilename(".")
    bamfile = P.snip(bamfile)
    statement = """
    samtools sort %(bamfile)s.bam -o %(T1)s.bam;
    samtools index %(T1)s.bam;
    mv %(T1)s.bam %(bamfile)s.bam;
    mv %(T1)s.bam.bai %(bamfile)s.bam.bai""" % locals()
    P.run()


def makeBamLink(currentname, newname):
    '''
    Make links to an existing bam file and its index
    '''
    cwd = os.getcwd()
    os.system("""
    ln -s %(cwd)s/%(currentname)s %(cwd)s/%(newname)s;
    ln -s %(cwd)s/%(currentname)s.bai %(cwd)s/%(newname)s.bai;
    """ % locals())

def makeLink(currentname, newname):
    '''
    Make links to an existing bam file and its index
    '''
    cwd = os.getcwd()
    os.system("""
    ln -s %(cwd)s/%(currentname)s %(cwd)s/%(newname)s;
    """ % locals())


def readTable(tabfile):
    '''
    Used to read the "peakcalling_bams_and_inputs.tsv" file back
    into memory.
    '''
    df = pd.read_csv(tabfile, sep="\t")
    chips = df['ChipBam'].values
    inputs = df['InputBam'].values
    pairs = zip(chips, inputs)
    D = dict(pairs)
    return D


class Peakcaller(object):
    ''' base class for peakcallers
    Peakcallers call peaks from a BAM file and then post-process peaks
    to generate a uniform output

    Attributes
    ----------

    threads: int
        number of threads to use for peakcalling
    paired_end: bool
        BAM contains paired end reads
    output_all: bool
        output all peaks (required for IDR)

    '''
    def __init__(self,
                 threads=1,
                 paired_end=True,
                 output_all=True,
                 tool_options=None):
        self.threads = threads
        self.paired_end = paired_end
        self.output_all = output_all
        self.tool_options = tool_options

    def callPeaks(self, infile, outfile, controlfile):
        ''' build command line statement to call peaks'''
        return ""

    def compressOutput(self, infile, outfile, contigsfile, controlfile):
        ''' build command line statement to compress outfiles'''
        return ""

    def postProcessPeaks(self, infile, outfile):
        ''' build command line statement to postprocess peaks'''
        return ""

    def preparePeaksForIDR(self, infile, outfile):
        ''' build command line statement to prepare the IDR input'''
        return ""

    def build(self, infile, outfile, contigsfile=None, controlfile=None,
              insertsizef=None, idr=0, idrc=0, idrsuffix=None, idrcol=None):
        ''' build complete command line statement'''

        peaks_outfile, peaks_cmd = self.callPeaks(infile, outfile, controlfile)
        compress_cmd = self.compressOutput(
            infile, outfile,  contigsfile, controlfile)
        postprocess_cmd = self.postProcessPeaks(
            infile, outfile, controlfile, insertsizef)

        if idr == 1:
            prepareIDR_cmd = self.preparePeaksForIDR(outfile, idrc, idrsuffix,
                                                     idrcol)
        else:
            prepareIDR_cmd = ""

        full_cmd = " checkpoint ;".join((
            peaks_cmd, compress_cmd, postprocess_cmd, prepareIDR_cmd))

        return full_cmd


class Macs2Peakcaller(Peakcaller):
    ''' call peaks with macs2
    '''

    def __init__(self,
                 threads=1,
                 paired_end=True,
                 output_all=True,
                 tool_options=None,
                 tagsize=None):
        super(Macs2Peakcaller, self).__init__(threads, paired_end,
                                              output_all, tool_options)
        self.tagsize = tagsize

    def callPeaks(self, infile,  outfile, controlfile=None):
        ''' build command line statement to call peaks.

        The _log file represents the logging output from MACS2

        The original location for the logging file (suffixed .macs2)
        is overwritten by the output from passing the macs2 xls to
        bed2table.py to give the final required outfile
        '''

        if self.tool_options:
            options = [self.tool_options]
        else:
            options = []
        if controlfile:
            options.append("--control %s" % controlfile)
        if self.tagsize is None:
            self.tagsize = BamTools.estimateTagSize(infile)

        options.append("--tsize %i" % self.tagsize)
        options = " ".join(options)

        # example statement: macs2 callpeak -t R1-paupar-R1.call.bam -c
        # R1-lacZ-R1.call.bam -f BAMPE -g 2.39e9 --verbose 5 --bw 150 -q
        # 0.01 -m 10 100000 --set-name test

        # format bam needs to be set explicitely, autodetection does not
        # work. Use paired end mode to detect tag size
        if self.paired_end:
            if not BamTools.isPaired(infile):
                raise ValueError(
                    "paired end has been specified but "
                    "BAM is not paired %" % infile)
            format_options = '--format=BAMPE'
        else:
            format_options = '--format=BAM'

        # --bdg --SPMR: ask macs to create a bed-graph file with
        # fragment pileup per million reads

        statement = '''
        macs2 callpeak
        %(format_options)s
        --treatment %(infile)s
        --verbose=10
        --name=%(outfile)s
        --qvalue=%%(macs2_max_qvalue)s
        --bdg
        --SPMR
        %(options)s
        >& %(outfile)s ;
        mv %(outfile)s %(outfile)s_log;
        ''' % locals()

        return outfile, statement

    def compressOutput(self, infile, outfile,
                       contigsfile, controlfile):
        ''' build command line statement to compress outfiles'''

        statement = []

        # compress macs bed files and index with tabix
        for suffix in ('peaks', 'summits'):
            bedfile = outfile + "_%s.bed" % suffix
            if os.path.exists(bedfile):
                statement.append('''
                     bgzip -f %(bedfile)s;
                     tabix -f -p bed %(bedfile)s.gz
                ''' % locals())

        # convert normalized bed graph to bigwig
        # saves 75% of space
        # compressing only saves 60%

        statement.append('''
        bedGraphToBigWig %(outfile)s_treat_pileup.bdg
        %(contigsfile)s %(outfile)s_treat_pileup.bw ;
        checkpoint ; rm -rf %(outfile)s_treat_pileup.bdg''' % locals())

        statement.append('''
        bedGraphToBigWig %(outfile)s_control_lambda.bdg
        %(contigsfile)s %(outfile)s_control_lambda.bw ;
        checkpoint ; rm -rf %(outfile)s_control_lambda.bdg''' % locals())

        # index and compress peak file
        suffix = 'peaks.xls'
        statement.append(
            '''grep -v "^$" < %(outfile)s_%(suffix)s
            | bgzip > %(outfile)s_%(suffix)s.gz;
            x=$(zgrep "[#|log]" %(outfile)s_%(suffix)s.gz | wc -l);
            tabix -f -b 2 -e 3 -S $x %(outfile)s_%(suffix)s.gz;
             checkpoint; rm -f %(outfile)s_%(suffix)s''' % locals())

        return "; checkpoint ;".join(statement)

    def postProcessPeaks(self, infile, outfile, controlfile, insertsizefile):
        '''
        postprocess MACS 2 results

        Output files from MACS2 are passed to bed2table to annotate
        them with various metrics including the centre of the peak and
        the number of reads in the sample bam and control bam at the
        peak position.

        The xls output from macs2 is passed to bed2table to generate
        the final outfile.

        Depending on the peak calling method of macs2,
        either the a summits or broadpeak macs2 output will also be passed to
        bed2table for annotation.
        '''

        filename_bed = outfile + "_peaks.xls.gz"
        filename_subpeaks = outfile + "_summits.bed"
        filename_broadpeaks = "%s_peaks.broadPeak" % outfile

        outfile_subpeaks = P.snip(
            outfile, ".macs2", ) + ".subpeaks.macs_peaks.bed"

        outfile_broadpeaks = P.snip(
            outfile, ".macs2", ) + ".broadpeaks.macs_peaks.bed"

        shift = getMacsPeakShiftEstimate(insertsizefile)
        assert shift is not None,\
            "could not determine peak shift from file %s" % insertsizefile

        peaks_headers = ",".join((
            "contig", "start", "end",
            "interval_id",
            "-log10\(pvalue\)", "fold", "-log10\(qvalue\)",
            "macs_nprobes", "macs_peakname"))

        if controlfile:
            control = '''--control-bam-file=%(controlfile)s
                         --control-offset=%(shift)s''' % locals()
        else:
            control = ""

        statement = '''
        zcat %(filename_bed)s |
        awk '$1!="chr"' |
        python %%(scriptsdir)s/bed2table.py
        --counter=peaks
        --bam-file=%(infile)s
        --offset=%(shift)i
        %(control)s
        --output-all-fields
        --output-bed-headers=%(peaks_headers)s
        --log=%(outfile)s.log
        > %(outfile)s ;
        ''' % locals()

        # check masc2 running mode
        if "--broad" in self.tool_options:

            broad_headers = ",".join((
                "contig", "start", "end",
                "interval_id",
                "Height"))

            # add a peak identifier and remove header
            statement += '''
            cat %(filename_broadpeaks)s |
            awk '/Chromosome/ {next; }
            {printf("%%%%s\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\n",
            $1,$2,$3,++a,$4)}'
            | python %%(scriptsdir)s/bed2table.py
            --counter=peaks
            --bam-file=%(infile)s
            --offset=%(shift)i
            %(control)s
            --output-all-fields
            --output-bed-headers=%(broad_headers)s
            --log=%(outfile)s.log
            > %(outfile_broadpeaks)s;
            ''' % locals()

        # if not broad, will generate subpeaks
        else:

            subpeaks_headers = ",".join((
                "contig", "start", "end",
                "interval_id",
                "Height"))

            # add a peak identifier and remove header
            statement += '''
            cat %(filename_subpeaks)s |
            awk '/Chromosome/ {next; }
            {printf("%%%%s\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\n",
            $1,$2,$3,++a,$5)}'
            | python %%(scriptsdir)s/bed2table.py
            --counter=peaks
            --bam-file=%(infile)s
            --offset=%(shift)i
            %(control)s
            --output-all-fields
            --output-bed-headers=%(subpeaks_headers)s
            --log=%(outfile)s.log
            > %(outfile_subpeaks)s ;''' % locals()

        return statement

    def preparePeaksForIDR(self, outfile, idrc, idrsuffix, idrcol):
        statement = ''
        idrout = "%s_IDRpeaks" % outfile
        narrowpeaks = "%s_peaks.%s" % (outfile, idrsuffix)
        col = idrcol
        statement += '''sort -h -r -k%(col)i,%(col)i %(narrowpeaks)s |
        head -%(idrc)s > %(idrout)s''' % locals()
        return statement


    def summarise(self, infile):
        '''Parses the MACS2 logfile to extract
        peak calling parameters and results.

        TODO: doesn't report peak numbers...

        Arguments
        ---------
        infiles : string
            Input files in MACS2 log file format.
        outfile : string
            Filename of output file in tab-delimited text format.
        '''

        infile = "%s_log" % infile
        outfile = "%s.table" % infile

        map_targets = [
            ("fragment size = (\d+)",
             "fragment_size", ()),
            ("total fragments in treatment:\s+(\d+)",
             "fragment_treatment_total", ()),
            ("fragments after filtering in treatment:\s+(\d+)",
             "fragment_treatment_filtered", ()),
            ("total fragments in control:\s+(\d+)",
             "fragment_control_total", ()),
            ("fragments after filtering in control:\s+(\d+)",
             "fragment_control_filtered", ())
        ]

        mapper, mapper_header = {}, {}
        for x, y, z in map_targets:
            mapper[y] = re.compile(x)
            mapper_header[y] = z

        keys = [x[1] for x in map_targets]

        results = collections.defaultdict(list)
        with IOTools.openFile(infile) as f:
            for line in f:
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".macs2_log")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                v = "\t".join(map(str, val + ["na"] * (c - len(val))))
            row.append(v)

        peaks = IOTools.openFile(
            infile.replace(".macs2_log",
                           ".macs2_peaks.xls.gz")).readlines()
        npeaks = 0
        for line in peaks:
            if "#" not in line and "log10" not in line:
                npeaks += 1

        row.extend([str(npeaks)])
        keys.extend(["number_of_peaks"])

        out = IOTools.openFile(outfile, "w")
        out.write("sample\t%s\n%s\n" % ("\t".join(keys), "\t".join(row)))
        out.close()


#############################################
# IDR Functions

@cluster_runnable
def makePairsForIDR(infiles, outfile, useoracle, df):
    pseudo_reps = []
    pseudo_pooled = []
    notpseudo_reps = []
    notpseudo_pooled = []

    for f in infiles:
        if "pseudo" in f and "pooled" in f:
            pseudo_pooled.append(f)
        elif "pseudo" in f:
            pseudo_reps.append(f)
        elif "pooled" in f:
            notpseudo_pooled.append(f)
        else:
            notpseudo_reps.append(f)

    oracledict = dict()
    cr_pairs = df['Condition'] + "_" + df['Tissue']
    for cr in cr_pairs:
        for npp in notpseudo_pooled:
            npp1 = npp.split("/")[-1]
            if npp1.startswith(cr):
                oracledict[cr] = npp
    i = 0
    for bam in df['bamReads']:
        bam = P.snip(bam)
        cr = cr_pairs[i]
        for npp in notpseudo_pooled:
            npp1 = npp.split("/")[-1]
            if npp1.startswith(cr):
                oracledict[bam] = npp
        i += 1

    pseudo_reps = np.array(sorted(list(set(pseudo_reps))))
    pseudo_pooled = np.array(sorted(list(set(pseudo_pooled))))

    into = len(pseudo_reps) / 2
    pseudoreppairs = np.split(pseudo_reps, into)
    pseudoreppairs = [tuple(item) for item in pseudoreppairs]
    pseudoreppairs_rows = []
    for tup in pseudoreppairs:
        stem = tup[0].split("/")[-1]
        stem = re.sub(r'_filtered.*', '', stem)
        oraclenam = oracledict[stem]
        if useoracle == 1:
            oracle = oraclenam
        else:
            oracle = "None"
        tissue, condition = oraclenam.split("/")[-1].split("_")[0:2]
        row = ((tup[0], tup[1], "self_consistency", oracle, tissue,
                condition))
        pseudoreppairs_rows.append(row)

    into = len(pseudo_pooled) / 2
    pseudopooledpairs = np.split(pseudo_pooled, into)
    pseudopooledpairs = [tuple(item) for item in pseudopooledpairs]
    pseudopooledpairs_rows = []
    for tup in pseudopooledpairs:
        stem = tup[0].split("/")[-1]
        stem = re.sub(r'_pooled_filtered.*', '', stem)
        oraclenam = oracledict[stem]
        if useoracle == 1:
            oracle = oraclenam
        else:
            oracle = "None"
        tissue, condition = oraclenam.split("/")[-1].split("_")[0:2]
        row = ((tup[0], tup[1], "pooled_consistency", oracle, tissue,
                condition))

        pseudopooledpairs_rows.append(row)

    # segregate the non pseudo non pooled files into sets of replicates
    # with the same condition and tissue
    conditions = df['Condition'].values
    tissues = df['Tissue'].values
    names = df['bamReads'].values
    i = 0
    repdict = dict()
    for c in conditions:
        t = tissues[i]
        n = names[i]
        repdict.setdefault("%s_%s" % (c, t), [])
        repdict["%s_%s" % (c, t)].append("%s_filtered" % P.snip(n))
        i += 1

    idrrepdict = dict()
    for k in repdict.keys():
        reps = repdict[k]
        idrrepdict.setdefault(k, [])
        for rep in reps:
            for f in notpseudo_reps:
                if rep in f:
                    idrrepdict[k].append(f)

    # generate every possible pair of replicates
    reppairs = []
    for k in idrrepdict.keys():
        reppairs += list(itertools.combinations(idrrepdict[k], 2))

    reppairs_rows = []
    for tup in reppairs:
        stem = tup[0].split("/")[-1]
        stem = re.sub(r'_filtered.*', '', stem)
        oraclenam = oracledict[stem]
        if useoracle == 1:
            oracle = oraclenam
        else:
            oracle = "None"
        tissue, condition = oraclenam.split("/")[-1].split("_")[0:2]
        row = ((tup[0], tup[1], "replicate_consistency", oracle, tissue,
                condition))
        reppairs_rows.append(row)
    pairs = pseudoreppairs_rows + pseudopooledpairs_rows + reppairs_rows

    out = IOTools.openFile(outfile, "w")
    out.write(
        "file1\tfile2\tIDR_comparison_type\toutput_file\tCondition\tTissue\n")
    for p in pairs:
        out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % p)

    out.close()


def buildIDRStatement(infile1, infile2, outfile,
                      sourcec, unsourcec, soft_idr_thresh, idrPARAMS, options,
                      oraclefile, test=False):
    statement = []
    statement.append(sourcec)

    log = "%s.log" % P.snip(outfile)

    inputfiletype = idrPARAMS['idrsuffix']
    rank = idrPARAMS['idrcolname']

    if oraclefile == "None":
        oraclestatement = ""
    else:
        oraclestatement = "--peak-list %s" % oraclefile

    if inputfiletype == "broadPeak":
        outputfiletype = "broadPeak"
    else:
        outputfiletype = "narrowPeak"

    T = P.getTempFilename(".")

    if test is True:
        idrstatement = """idr --version >>%(log)s;
                          idr
                          --samples %(infile1)s %(infile2)s
                          --output-file %(outfile)s
                          --soft-idr-threshold %(soft_idr_thresh)s
                          --input-file-type %(inputfiletype)s
                          --rank %(rank)s
                          --output-file-type %(outputfiletype)s
                           %(oraclestatement)s
                           %(options)s
                          --only-merge-peaks
                          --verbose 2>>%(log)s
        """ % locals()

    else:
        idrstatement = """idr --version >>%(log)s;
                          idr
                          --samples %(infile1)s %(infile2)s
                          --output-file %(outfile)s
                          --plot
                          --soft-idr-threshold %(soft_idr_thresh)s
                          --input-file-type %(inputfiletype)s
                          --rank %(rank)s
                          --output-file-type %(outputfiletype)s
                           %(oraclestatement)s
                           %(options)s
                          --verbose 2>>%(log)s
        """ % locals()
    statement.append(idrstatement)
    statement.append(unsourcec)
    statement = "; ".join(statement)
    return statement


def summariseIDR(infiles, outfile, pooledc, selfc, repc):
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
                                                  'Tissue',
                                                  'IDR_Successful'])

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
            thresholds.append(pooledc)
        elif item == "self_consistency":
            thresholds.append(selfc)
        elif item == "replicate_consistency":
            thresholds.append(repc)

    alltab['IDR_Thresholds'] = thresholds
    alltab['IDR_Thresholds'] = alltab['IDR_Thresholds'].astype('float')
    alltab['IDR_Transformed_Thresholds'] = -(np.log10(thresholds))

    alltab = alltab.rename(columns={'IDR_Successful_x': 'IDR_Successful'})
    alltab = alltab[['Output_Filename', 'Replicate_Type', 'Total_Peaks',
                     'Peaks_Passing_IDR', 'Percentage_Peaks_Passing_IDR',
                     'Peaks_Failing_IDR', 'IDR_Thresholds',
                     'IDR_Transformed_Thresholds',
                     'Input_File_1', 'Input_File_2', 'Condition', 'Tissue',
                     'Oracle_Peak_File', 'IDR_Successful']]

    alltab['Condition'] = alltab['Condition'].astype('str')
    alltab['Tissue'] = alltab['Tissue'].astype('str')
    alltab['Experiment'] = alltab['Condition'] + "_" + alltab['Tissue']
    alltab['Conservative_Peak_List'] = False
    replicates = alltab[alltab['Replicate_Type'] == 'replicate_consistency']
    for exp in set(replicates['Experiment'].values):
        subtab = replicates[replicates['Experiment'] == exp]
        m = max(subtab['Peaks_Passing_IDR'])
        val = subtab[subtab['Peaks_Passing_IDR'] == m].index.values
        alltab.ix[val, 'Conservative_Peak_List'] = True

    alltab['Optimal_Peak_List'] = False
    for exp in set(alltab['Experiment'].values):
        subtab = alltab[((alltab['Experiment'] == exp) &
                         ((alltab['Replicate_Type'] == 'pooled_consistency') |
                          (alltab['Replicate_Type'] ==
                           'replicate_consistency')))]
        m = max(subtab['Peaks_Passing_IDR'])
        val = subtab[subtab['Peaks_Passing_IDR'] == m].index.values
        alltab.ix[val, 'Optimal_Peak_List'] = True

    alltab.to_csv(outfile, sep="\t", index=False)


def doIDRQC(infile, outfile):
    alltab = pd.read_csv(infile, sep="\t")
    outputrows = []
    self_c_alltab = alltab[alltab['Replicate_Type'] == 'self_consistency']
    pooled_c_alltab = alltab[alltab['Replicate_Type'] == "pooled_consistency"]
    replicate_c_alltab = alltab[alltab['Replicate_Type'] ==
                                'replicate_consistency']

    for exp in set(self_c_alltab['Experiment'].values):

        # Self-Consistency
        # Make a sub-table of only the self consistency IDRs for this
        # experiment
        # (Tissue and Condition)
        subtab = self_c_alltab[self_c_alltab['Experiment'] == exp]

        # Take the number of peaks passing IDR for each replicate of the
        # experiment
        # and the names of the files containing these peaks
        Ns = subtab['Peaks_Passing_IDR'].values
        reps = subtab['Output_Filename'].values

        # Generate every possible pair of replicates within the experiment
        pairs = list(itertools.combinations(range(len(reps)), 2))

        names = []

        # For each pair of replicates
        for pair in pairs:
            # record the filenames of this pair
            names = ((reps[pair[0]], reps[pair[1]]))

            # record the number of peaks passing IDR for this pair
            N1N2 = ((Ns[pair[0]], Ns[pair[1]]))

            m = max(N1N2)
            i = N1N2.index(m)

            maxN1N2 = float(N1N2[i])
            maxname = names[i]
            pmax = pair[i]

            if i == 0:
                minN1N2 = float(N1N2[1])
                minname = names[1]
                pmin = pair[1]
            else:
                minN1N2 = float(N1N2[0])
                minname = names[0]
                pmin = pair[0]

            # calculate the self-consistency ratio -
            # the maxiumum of N1 and N2 divided by the minimum of N1 and N2
            # This should be less than 2
            if minN1N2 == 0:
                self_con_ratio = float('nan')
            else:
                self_con_ratio = maxN1N2 / minN1N2

            # store the self-consistency ratio for this pair
            subname1 = re.sub("_filtered.*", "", maxname)
            subname2 = re.sub("_filtered.*", "", minname)

            row = ((exp, 'self-consistency ratio (N%i / N%i)' % (pmax, pmin),
                    '%s / %s' % (subname1, subname2), self_con_ratio,
                    names[0], names[1]))

            outputrows.append(row)

        # Comparison Pooled Pseudo-replicates

        # Take the subset of the table containing pooled pseudoreplicates
        # for this
        # experiment
        pooledsubtab = pooled_c_alltab[pooled_c_alltab['Experiment'] == exp]

        # Find Np - the number of peaks which pass the IDR threshold when
        # comparing
        # pooled pseudoreplicates
        Np = pooledsubtab['Peaks_Passing_IDR'].values[0]

        # Keep track of which file this is
        poolednam = pooledsubtab['Output_Filename'].values[0]
        pooledsubname = re.sub("_filtered.*", "", poolednam)

        # Take the subset of the table containing true replicates for this
        # experiment
        replicatesubtab = replicate_c_alltab[
            replicate_c_alltab['Experiment'] == exp]

        # Extract the names and the number of peaks passing IDR from this
        # table
        replicatenames = list(replicatesubtab['Output_Filename'].values)
        replicatescores = list(replicatesubtab['Peaks_Passing_IDR'].values)

        # Find Nt - the maximum number of peaks passing IDR for true replicates
        Nt = max(replicatescores)
        whichNt = replicatescores.index(Nt)

        # Keep track of which file Nt came from
        maxrepname = replicatenames[whichNt]
        maxrepsubname = re.sub("_filtered.*", "", maxrepname)

        # The rescue ratio is max(Np, Nt) / min (Np, Nt) - find which
        # is max and which is min
        NpNt = ((Np, Nt))
        maxNpNt = max(NpNt)
        minNpNt = min(NpNt)

        # Rescue ratio
        if min(NpNt) == 0:
            RR = float('nan')
        else:
            RR = float(maxNpNt) / float(minNpNt)

        i = NpNt.index(maxNpNt)
        if i == 0:
            NpNtstring = "Np / Nt"
            f1 = pooledsubname
            f2 = maxrepsubname
        else:
            NpNtstring = "Nt / Np"
            f1 = maxrepsubname
            f2 = pooledsubname

        formula = "%s / %s" % (f1, f2)

        # Generate a table row with this data
        row = ((
            exp, 'rescue ratio (%s)' % NpNtstring, formula, RR, ",".join(
                replicatenames), poolednam))

        # keep this data to put in the output table
        outputrows.append(row)

    outputtab = pd.DataFrame(outputrows,
                             columns=['Experiment', 'Ratio_Type',
                                      'Ratio_Formula',  'Ratio',
                                      'Comparison_File_1',
                                      'Comparison_File_2'])

    outputtab[
        'Individual_Reproducibility'] = [
            'PASS' if x < 2 else 'WARN' for x in outputtab['Ratio']]

    outputtab['Experiment_Reproducibility'] = ""

    finaltab = pd.DataFrame(columns=outputtab.columns)

    for exp in set(outputtab['Experiment'].values):
        exptab = outputtab[outputtab['Experiment'] == exp]
        reps = list(set(exptab['Individual_Reproducibility']))
        if len(reps) != 1:
            exptab['Experiment_Reproducibility'] = "WARN"
        elif reps[0] == "PASS":
            exptab['Experiment_Reproducibility'] = "PASS"
        elif reps[0] == "WARN":
            exptab['Experiment_Reproducibility'] = "FAIL"
        finaltab = finaltab.append(exptab)

    finaltab.to_csv(outfile, sep="\t", index=None)

#############################################
# QC Functions


def runCHIPQC(infile, outfiles, rdir):
    runCHIPQC_R = R('''
    function(samples, outdir, cwd){
        library("ChIPQC")
        cwd = "/ifs/projects/katherineb/test_data_peakcalling/test3"
        print ("cwd")
        print (cwd)
        print ("samples")
        print (samples)
        setwd(cwd)
        print (getwd())
        samples$Tissue = as.factor(samples$Tissue)
        samples$Factor = as.factor(samples$Factor)
        samples$Replicate = as.factor(samples$Replicate)
        experiment = ChIPQC(samples)
        ChIPQCreport(experiment, reportFolder=outdir)
    }
    ''' % locals())
    cwd = os.getcwd()
    runCHIPQC_R(pandas2ri.py2ri(infile), rdir, cwd)


def summariseMACS2(infiles, outfile):
    '''Parses the MACS2 logfile to extract
    peak calling parameters and results.

    TODO: doesn't report peak numbers...

    Arguments
    ---------
    infiles : string
        Input files in MACS2 log file format.
    outfile : string
        Filename of output file in tab-delimited text format.
    '''

    def __get(line, stmt):
        x = line.search(stmt)
        if x:
            return x.groups()

    # mapping patternts to values.
    # tuples of pattern, label, subgroups
    map_targets = [
        ("tags after filtering in treatment:\s+(\d+)",
         "fragment_treatment_filtered", ()),
        ("total tags in treatment:\s+(\d+)",
         "fragment_treatment_total", ()),
        ("tags after filtering in control:\s+(\d+)",
         "fragment_control_filtered", ()),
        ("total tags in control:\s+(\d+)",
         "fragment_control_total", ()),
        ("predicted fragment length is (\d+) bps",
         "fragment_length", ()),
        # Number of peaks doesn't appear to be reported!.
        ("#3 Total number of candidates: (\d+)",
         "ncandidates", ("positive", "negative")),
        ("#3 Finally, (\d+) peaks are called!",
         "called",
         ("positive", "negative"))
    ]

    mapper, mapper_header = {}, {}
    for x, y, z in map_targets:
        mapper[y] = re.compile(x)
        mapper_header[y] = z

    keys = [x[1] for x in map_targets]

    outs = IOTools.openFile(outfile, "w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend(["%s_%s" % (k, x) for x in mapper_header[k]])
        else:
            headers.append(k)
    outs.write("track\t%s" % "\t".join(headers) + "\n")

    for infile in infiles:
        results = collections.defaultdict(list)
        with IOTools.openFile(infile) as f:
            for line in f:
                if "diag:" in line:
                    break
                for x, y in mapper.items():
                    s = y.search(line)
                    if s:
                        results[x].append(s.groups()[0])
                        break

        row = [P.snip(os.path.basename(infile), ".macs2")]
        for key in keys:
            val = results[key]
            if len(val) == 0:
                v = "na"
            else:
                c = len(mapper_header[key])
                # append missing data (no negative peaks without control files)
                v = "\t".join(map(str, val + ["na"] * (c - len(val))))
            row.append(v)
            # assert len(row) -1 == len( headers )
        outs.write("\t".join(row) + "\n")

    outs.close()
