'''
======================================================================

Requirements:


Reference
---------

'''

import os
import re
import collections
import pandas
import math
import numpy
import numpy.ma as ma
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

##############################################
# Preprocessing Functions


def trackFilters(filtername, bamfile, tabout):
    return """echo %(filtername)s >> %(tabout)s;
              samtools view -c %(bamfile)s.bam >> %(tabout)s; """ % locals()


def appendSamtoolsFilters(statement, inT, tabout, filters, pe):
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
        elif filt == "unpaired":
            if pe is True:
                string = "-f 1 -f 2"
            else:
                E.warn("""Bam file is not paired-end so unpaired reads have
                not been filtered""")
                continue
        # append any new samtools filters to this list
        if filt in ["unmapped", "unpaired"]:
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


def filterBams(infile, outfiles, filters, bedfiles, blthresh, pe, strip,
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
    statement, inT = appendSamtoolsFilters(statement, inT, tabout, filters, pe)
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
def makePseudoBams(infile, outfiles, pe, randomseed):
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
            if j % 2 == 0:
                destint = intlist[i]
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
    #   pysam.sort(outf, T, catch_stdout=False)
    #   pysam.index("%s.bam" % T, catch_stdout=False)
        bamfile = pysam.AlignmentFile("%s.bam" % T, "rb")
        all = []
        for nam in bamfile:
            all.append(nam.qname)
        if pe is True:
            assert (
                len(set(all)) ==
                len(all) / 2), "Error splitting bam file %(outf)s" % locals()
        lens.append(bamfile.count())
    #     os.remove(T)
        shutil.move("%s.bam" % T, outf)
        shutil.move("%s.bam.bai" % T, "%s.bai" % outf)
    E.info("Bamfile 1 length %i, Bamfile 2 length %i" % (lens[0], lens[1]))

#############################################
# Peakcalling Functions


# TS - move to top of file later:
def getMacsPeakShiftEstimate(infile):
    ''' parse the peak shift estimate file from macs,
    return the fragment size '''

    with IOTools.openFile(infile, "r") as inf:

        header = inf.next().strip().split("\t")
        values = inf.next().strip().split("\t")

        fragment_size_mean_ix = header.index("fragmentsize_mean")

        fragment_size = int(float(values[fragment_size_mean_ix]))

        return fragment_size


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

    def build(self, infile, outfile, contigsfile=None, controlfile=None,
              insertsizef=None):
        ''' build complete command line statement'''

        peaks_outfile, peaks_cmd = self.callPeaks(infile, outfile, controlfile)
        compress_cmd = self.compressOutput(
            infile, outfile,  contigsfile, controlfile)
        postprocess_cmd = self.postProcessPeaks(
            infile, outfile, controlfile, insertsizef)

        full_cmd = " checkpoint ;".join((
            peaks_cmd, compress_cmd, postprocess_cmd))

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
        the _raw file represents the unadulterated output from
        MACS2. The actual output is then processed by the
        postProcess method.  This _raw file is retained for
        reproducibility.
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
        mv %(outfile)s %(outfile)s_raw;
        ''' % locals()

        return outfile, statement

    def compressOutput(self, infile, outfile,
                       contigsfile, controlfile):
        ''' build command line statement to compress outfiles

        ### NOTE
        This takes a contig file but I don't understand what it is doing
        with it.
        The contig file was hard coded - I have moved it to the pipeline.ini.
        '''

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
            tabix -f -b 2 -e 3 -S 26 %(outfile)s_%(suffix)s.gz;
             checkpoint; rm -f %(outfile)s_%(suffix)s''' % locals())

        return "; checkpoint ;".join(statement)

    def postProcessPeaks(self, infile, outfile, controlfile, insertsizefile):
        '''
        postprocess MACS 2 results


        ### NOTE
        This is doing a subset of the steps the old peakcalling pipeline
        did on the MACS2 output - we're not sure how you chose which are still
        needed and which to delete?

        Some of the steps look slightly dodgy - the step which is generating
        the subpeaks seems to just take the "summits" bed file and write
        it out again but convert all the columns, including those which are
        strings, into integers?
        ###

        '''

        filename_bed = outfile + "_peaks.xls.gz"
        filename_diag = outfile + "_diag.xls"
        filename_r = outfile + "_model.r"
        filename_rlog = outfile + ".r.log"
        filename_pdf = outfile + "_model.pdf"
        filename_subpeaks = outfile + "_summits.bed"

        outfile_subpeaks = P.snip(
            outfile, ".macs2", ) + ".subpeaks.macs_peaks.bed"

        filename_broadpeaks = P.snip(outfile,
                                     ".macs2") + ".broadpeaks.macs_peaks.bed"
        outfile_broadpeaks = P.snip(
            outfile, ".macs2", ) + ".broadpeaks.macs_peaks.bed"

        # TS - will bedfile ever not be created?
        # previous pipeline had a test to check. IS there an option
        # which prevents bedfile output?
        # PREVIOUS CODE:
        # if not os.path.exists(filename_bed):
        #    E.warn("could not find %s" % filename_bed)
        #    P.touch(outfile)
        #    return

        shift = getMacsPeakShiftEstimate(insertsizefile)
        assert shift is not None,\
            "could not determine peak shift from file %s" % insertsizefile

        peaks_headers = ",".join((
            "contig", "start", "end",
            "interval_id",
            "pvalue", "fold", "qvalue",
            "macs_nprobes"))

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
                "Height",
                "SummitPosition"))

            # add a peak identifier and remove header
            statement += '''
            cat %(filename_subpeaks)s |
            awk '/Chromosome/ {next; }
            {printf("%%%%s\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\n",
            $1,$2,$3,++a,$4,$5)}'
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


#############################################
# QC Functions
