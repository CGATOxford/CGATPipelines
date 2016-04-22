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

##############################################
# Preprocessing Functions

def filterBams(infile, outfiles, filters, pe):
    bamout, out = outfiles
    tabout = open(out, "w")
    tabout.write("unmapped\t%s\n" % ("\t".join(filters)))
    tabout.close()

    F = []
    if 'unmapped' in filters:
        F.append(str(4))
    if 'unpaired' in filters and pe is True:
        F.append(str(12))
    else:
        E.warn("""Bam file is not paired-end so unpaired reads have not been
        filtered""")

    inT = infile
    outT = P.getTempFilename(".")

    # ensure bamfile is sorted,
    statement = """samtools sort %(inT)s %(outT)s; """ % locals()
    statement += """samtools view -c %(outT)s.bam >> %(out)s; """ % locals()


    tempfiles = []
    i = 0

    for num in F:
        T = P.getTempFilename(".")
        if i == 0:
            inT = outT
            outT = P.getTempFilename(".")
        else:
            inT = outT
            outT = P.getTempFilename(".")

        statement += """samtools view -b -F %(num)s %(inT)s.bam
        > %(outT)s.bam; """ % locals()
        statement += """samtools view -c %(outT)s.bam >> %(out)s; """ % locals()
        tempfiles.append(T)
        i += 1

    # remove duplicates, if requested
    inT = outT
    outT = P.getTempFilename(".")
    if 'duplicates' in filters:
        statement += """
        MarkDuplicates \
        INPUT=%(inT)s.bam \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        OUTPUT=%(outT)s.bam \
        METRICS_FILE=/dev/null \
        VALIDATION_STRINGENCY=SILENT \
        2> %(outT)s.log; """ % locals()

        statement += """samtools view -c %(outT)s.bam >> %(out)s; """ % locals()

    statement = statement.replace("\n", "")
    print statement
#    P.run()


def estimateInsertSize(infiles, outfile, insert_params):
    res = []
    for infile in infiles:
        fragment_length = BamTools.estimateInsertSizeDistribution(
            infile,
            alignments=insert_params['alignments'],
            n=insert_params['n'],
            method='picard',
            similarity_threshold=insert_params['st'])
        res.append(fragment_length)
    res = pd.DataFrame(res, columns=['mean_insert_size', 'standard_dev',
                                     'pairs_used_for_estimation'])
    res.to_csv(outfile, sep="\t")

#############################################
# Peakcalling Functions


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

        peaks_cmd = self.callPeaks(infile, outfile, controlfile)
        compress_cmd = self.compressOutput(
            infile, outfile, contigsfile, controlfile)
        postprocess_cmd = self.postProcessPeaks(infile, outfile)

        full_cmd = "; checkpoint ;".join((peaks_cmd, compress_cmd, postprocess_cmd))

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
        >& %%(outfile)s
        ''' % locals()

        return statement

    def compressOutput(self, infile, outfile, contigsfile, controlfile=None):
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
            tabix -f -b 2 -e 3 -S 26 %(outfile)s_%(suffix)s.gz;
             checkpoint; rm -f %(outfile)s_%(suffix)s''' % locals())

        return "; checkpoint ;".join(statement)




#############################################
# QC Functions
