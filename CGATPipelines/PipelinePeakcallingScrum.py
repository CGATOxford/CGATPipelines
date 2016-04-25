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

##############################################
# Preprocessing Functions
def filterBams(infile, outfile, filters):

    F = ""
    if 'unmapped' in filters:
        F += str(4)
    if 'unpaired' in filters:
        F += str(12)

    T1 = P.getTempFilename(".")
    T2 = P.getTempFilename(".")
    # ensure bamfile is sorted,
    statement = """samtools sort %(infile)s %(T1)s;
           samtools index %(T1)s;
           samtools view %(T1)s -b -F %(F)s > %(T2)s;"""

    # remove duplicates, if requested
    if 'duplicates' in filters:
        statement += """
        MarkDuplicates \
        INPUT=%(T2)s \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        OUTPUT=%(outfile)s \
        METRICS_FILE=/dev/null \
        VALIDATION_STRINGENCY=SILENT \
        2> %(outfile)s.log"""

    P.run()

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

    def build(self, infile, outfile, contigsfile=None, controlfile=None):
        ''' build complete command line statement'''

        peaks_outfile, peaks_cmd = self.callPeaks(infile, outfile, controlfile)
        compress_cmd = self.compressOutput(
            infile, outfile,  contigsfile, controlfile)
        postprocess_cmd = self.postProcessPeaks(
            infile, outfile, controlfile)

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
        ''' build command line statement to call peaks'''

        if self.tool_options:
            options = [self.tool_options]
        else:
            options = []
        if controlfile:
            options.append("--control %s" % controlfile)
        if self.tagsize is not None:
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
                    "paired end has been specified but BAM is not paired %" % infile)
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
        ''' build command line statement to compress outfiles'''

        statement = []

        # compress macs bed files and index with tabix
        for suffix in ('peaks', 'summits'):
            bedfile = outfile + "_%s.bed" % suffix
            if os.path.exists(bedfile):
                statement.append('''
                     bgzip -f %(bedfile)s;
                     tabix -f -p bed %(bedfile)s.gz;
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

    def postProcessPeaks(self, infile, outfile, controlfile):
        ''' build command line statement to postprocess peaks'''

        filename_bed = outfile + "_peaks.xls.gz"
        filename_diag = outfile + "_diag.xls"
        filename_r = outfile + "_model.r"
        filename_rlog = outfile + ".r.log"
        filename_pdf = outfile + "_model.pdf"

        filename_subpeaks = outfile + "_summits.bed"
        outfile_subpeaks = P.snip(
            outfile, ".macs2", ) + ".subpeaks.macs_peaks.bed"

        filename_broadpeaks = P.snip(outfile, ".macs2") + ".broadpeaks.macs_peaks.bed"
        outfile_broadpeaks = P.snip(
            outfile, ".macs2", ) + ".broadpeaks.macs_peaks.bed"

        # TS - will bedfile ever not be created?
        # previous pipeline had a test to check. IS there an option
        # which prevents bedfile output?
        # PREVIOUS CODE:
        #if not os.path.exists(filename_bed):
        #    E.warn("could not find %s" % filename_bed)
        #    P.touch(outfile)
        #    return

        # TS - hardcodes "fragment.dir" and ".fragment_size" suffix
        # another option is to set fragment_dir and fragment suffix attributes

        peak_est_inf = os.path.join(
            "fragment_size.dir", "%s.fragment_size" % P.snip(infile, ".genome.bam"))
        shift = getMacsPeakShiftEstimate(peak_est_inf)
        assert shift is not None,\
            "could not determine peak shift from file %s" % peak_est_inf

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
            awk '/Chromosome/ {next; } {printf("%%%%s\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\n", $1,$2,$3,++a,$4)}'
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
            awk '/Chromosome/ {next; } {printf("%%%%s\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\t%%%%i\\n", $1,$2,$3,++a,$4,$5)}'
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
