'''
PipelinePreprocess.py - Utility functions for processing short reads
====================================================================

:Author: Tom smith
:Release: $Id$
:Date: |today|
:Tags: Python

UPDATE: Processing reads is a common task before mapping. Different sequencing
technologies and applications require different read processing.
This module provides utility functions to abstract some of these variations

The pipeline does not know what kind of data it gets (a directory
might contain single end or paired end data or both).

The module currently provides modules to perform:
    * hard-trimming (fastx_trimmer, trimmomatic)
    * adapter trimming (trimgalore, trimmomatic, cutadapt)
    * read end quality trimming (trimgalore, trimmomatic)
    * sliding window quality trimming (sickle, trimmomatic)
    * RRBS-specific trimming (trimgalore)

It has been tested with:
   * .fastq: paired-end and single-end



To do
=====

How to deal with flash?
e.g this is the only tool that doesn't do one-to-one processing

Code
----

'''

import re
import os
import shutil
import collections
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq
import CGAT.Genomics as Genomics
import CGATPipelines.PipelineMapping as Mapping
import pandas.io.sql as pdsql
import CGAT.Sra as Sra

SequenceInformation = collections.namedtuple("SequenceInformation",
                                             """paired_end
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""")


def getReadLengthFromFastq(filename):
    '''return readlength from a fasta/fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.

    '''

    with IOTools.openFile(filename) as infile:
        record = iterate(infile).next()
        readlength = len(record.seq)
        return readlength


def reverseComplement(sequence):
    '''
    Reverse complement a sequence
    '''
    sequence = sequence.upper()
    rev_seq = ''
    for nt in sequence:
        if nt == 'A':
            rev_seq += 'T'
        elif nt == 'T':
            rev_seq += 'A'
        elif nt == 'C':
            rev_seq += 'G'
        elif nt == 'G':
            rev_seq += 'C'
        elif nt == 'N':
            rev_seq = 'N'

    rev_seq = rev_seq[::-1]

    return rev_seq


def makeAdaptorFasta(infile, dbh, contaminants_file, outfile):
    '''
    Generate a .fasta file of adaptor sequences that are overrepresented
    in the reads from a sample.
    Requires cutadapt >= 1.7
    '''

    sample = infile.split("/")[-1].rstrip(".gz")
    # replace '-'  and '.' with '_'
    sample = sample.replace("-", "_")
    sample = sample.replace(".", "_")
    sample = sample.split("_")

    # handle fastq, paired-end and sra
    if infile.endswith(".fastq.gz"):
        sample.remove("fastq")
    elif infile.endswith(".sra"):
        sample.remove("sra")
        # depending on whether sra contains single or paired-end data
        # different behaviour is implemented
        outdir = P.getTempDir()
        f, format = Sra.peek(infile, outdir)
        E.info("sra file contains the following files: %s" % f)
        shutil.rmtree(outdir)
        if len(f) == 1:
            pass
        elif len(f) == 2:
            sample2 = sample+["2"]
            sample += ["1"]
    elif infile.endswith(".fastq.1.gz"):
        # sample.remove("fastq")
        pass
    elif infile.endswith(".fastq.2.gz"):
        # sample.remove("fastq")
        pass
    else:
        raise AttributeError("unrecognised sequence file format")

    sample = "_".join(sample)
    query = "SELECT * FROM %s_fastqc_Overrepresented_sequences;" % sample
    df = pdsql.read_sql(query, dbh, index_col=None)

    # this section handles paired end data from sra files
    if len(f) == 2:
        sample2 = "_".join(sample2)
        query = "SELECT * FROM %s_fastqc_Overrepresented_sequences;" % sample2
        df2 = pdsql.read_sql(query, dbh, index_col=None)
        df = df.append(df2)
        df = df.reset_index(drop=True)

    # if there are no over represented sequences break here
    if not len(df):
        P.touch(outfile)
        return None

    # generate contamination sequence dictionary
    contam_dict = {}
    for idx in df.index:
        overid = df.loc[idx]['Possible_Source']
        overid = overid.split(" (")[0]
        if not re.search("No Hit", overid):
            overid = overid.replace(",", "")
            overid = overid.replace(" ", "_")
            seqid = df.loc[idx]['Sequence']
            contam_dict[overid] = seqid
        else:
            pass

    overreps = set(df['Possible_Source'])
    # pull out suspected adaptor contamination name
    ids = [x.split(" (")[0] for x in overreps if not re.search("No Hit", x)]
    ids = [h.replace(",", "") for h in ids]
    ids = [g.replace(" ", "_") for g in ids]

    # line comments begin with '#'
    # adaptor name and sequence are split with tabs and end with both newline
    # and carriage returns.  Put these in a reference dictionary

    with IOTools.openFile(contaminants_file, "r") as cfile:
        lines = cfile.readlines()
    adapt_list = [l for l in lines if not re.search("#", l)]
    adapt_dict = {}

    for each in adapt_list:
        # source of bugs - row names contains whitespace
        # that may interfere with down stream processing of fasta file
        # remove extraneous ','
        each = each.split("\t")
        each = [k.rstrip("\r\n") for k in each]
        each = [h for h in each if len(h)]
        if len(each):
            seq_id = each[0].replace(",", "")
            seq_id = seq_id.replace(" ", "_")
            adapt_dict[seq_id] = each[1]
        else:
            pass

    with IOTools.openFile(outfile, "w") as ofile:
        for ad in ids:
            try:
                ofile.write(">%s\n%s\n" % (ad, adapt_dict[ad]))
            except KeyError:
                ofile.write(">%s\n%s\n" % (ad, contam_dict[ad]))


def mergeAdaptorFasta(infiles, outfile):
    '''
    Merge fasta files of adapter contamination,
    include reverse complement, remove duplicate sequences
    '''

    fasta_dict = {}
    for each in infiles:
        with IOTools.openFile(each, "r") as infle:
            for line in infle:
                if line[0] == '>':
                    adapt = line.lstrip(">").rstrip("\n")
                    fasta_dict[adapt] = set()
                    fasta_dict[adapt + "_R"] = set()
                else:
                    seq = line.rstrip("\n")
                    rev_seq = Genomics.complement(seq)
                    fasta_dict[adapt].add(seq)
                    fasta_dict[adapt + "_R"].add(rev_seq)

    # if there are no adapters to remove break the pipeline here
    if not len(fasta_dict):
        raise AttributeError("There are no overrepresented sequences in "
                             "these fastq files.  Please turn off this "
                             "feature and re-run the pipeline")
    else:
        pass
    with IOTools.openFile(outfile, "w") as outfle:
        for key, value in fasta_dict.items():
            outfle.write(">%s\n%s\n" % (key, list(value)[0]))


def getSequencingInformation(track):
    '''glean sequencing information from *track*.'''

    colour = False
    if os.path.exists("%s.fastq.gz" % track):
        first_pair = "%s.fastq.gz" % track
        second_pair = None
    elif os.path.exists("%s.fastq.1.gz" % track):
        first_pair = "%s.fastq.1.gz" % track
        second_pair = "%s.fastq.2.gz" % track
    elif os.path.exists("%s.csfasta.gz" % track):
        first_pair = "%s.csfasta.gz" % track
        second_pair = None
        colour = True

    second_length = None
    if second_pair:
        if not os.path.exists(second_pair):
            raise IOError("could not find second pair %s for %s" %
                          (second_pair, first_pair))
        second_length = getReadLengthFromFastq(second_pair)

    return SequenceInformation._make((second_pair is not None,
                                      first_pair, second_pair,
                                      getReadLengthFromFastq(first_pair),
                                      second_length,
                                      colour))


class MasterProcessor(Mapping.SequenceCollectionProcessor):

    '''Processes reads with tools specified.
    All in a single statement to be send to the cluster.
    '''

    # compress temporary fastq files with gzip
    compress = False

    def __init__(self,
                 save=True,
                 summarize=False,
                 threads=1,
                 *args, **kwargs):
        self.save = save
        self.summarize = summarize
        self.threads = threads
        if self.save:
            self.outdir = "processed.dir"
        else:
            self.outdir = P.getTempDir("/ifs/scratch")

        self.processors = []

    def add(self, processor):
        """add a processor to the list of tools to be executed."""
        self.processors.append(processor)

    def cleanup(self):
        '''clean up.'''
        if self.save:
            statement = 'checkpoint;'
        else:
            statement = 'rm -rf %s;' % self.outdir
        statement += '''rm -rf %s;''' % (self.tmpdir_fastq)
        return statement

    def build(self, infile, output_prefix, track):
        '''build a preprocessing statement

        Arguments
        ---------
        track : string
            prefix for output files. The suffix ".fastq.gz" will
            be appended for single-end data sets. For paired end
            data sets, two files will be created ending in
            ".fastq.1.gz" and ".fastq.2.gz"

        '''

        # preprocess to convert non-fastq into fastq
        # if already fastq, outfiles will be infiles
        cmd_preprocess, current_files = self.preprocess(infile, track + ".gz")

        # preprocess currently returns a list of lists
        current_files = current_files[0]

        cmd_processors = []

        # build statements from processors
        for idx, processor in enumerate(self.processors):

            # number of output files created in this step
            nfiles = processor.get_num_files(current_files)

            if nfiles == 2:
                suffixes = [track + ".fastq.1.gz",
                            track + ".fastq.2.gz"]
            else:
                suffixes = [track + ".fastq.gz"]

            if idx == len(self.processors)-1:
                # last iteration, write to output files
                next_files = [output_prefix + s
                              for s in suffixes]
            else:
                # add a prefix to each file denoting the processor
                next_files = [processor.prefix + "-" + track + s
                              for s in suffixes]

                # add scratch temporary directory to next files
                next_files = [os.path.join(self.tmpdir_fastq, x)
                              for x in next_files]

            cmd_process = processor.build(
                current_files,
                next_files,
                output_prefix + track + "-" + processor.prefix)

            current_files = next_files
            cmd_processors.append(cmd_process)

            if self.summarize:
                for fn in current_files:
                    cmd_processors.append(
                        """zcat %(fn)s
                        | python %%(scriptsdir)s/fastq2summary.py
                        --guess-format=illumina-1.8
                        > %(fn)s.summary""")

        cmd_process = " checkpoint; ".join(cmd_processors)
        cmd_clean = self.cleanup()

        assert cmd_preprocess.strip().endswith(";")
        assert cmd_process.strip().endswith(";")
        assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_process,
                                          cmd_clean))
        return statement


class ProcessTool(object):
    '''defines class attributes for a sequence utility tool'''

    def __init__(self, options, threads=1):
        self.processing_options = options
        self.threads = threads

    def get_num_files(self, infiles):
        """return the number of outputfiles created by
        processing `infiles`"""

        # the default is that the same number of files are
        # returned
        return len(infiles)


class Trimgalore(ProcessTool):

    prefix = 'trimgalore'

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == len(outfiles)
        assert len(infiles) in (1, 2)

        offset = Fastq.getOffset("sanger", raises=False)
        processing_options = self.processing_options

        if len(infiles) == 1:
            infile = infiles[0]
            outfile = outfiles[0]
            trim_out = "%s_trimmed.fq.gz" % (output_prefix)
            cmd = '''trim_galore %(processing_options)s
            --phred%(offset)s
            --output_dir %(outdir)s
            %(infile)s
            2>>%(output_prefix)s.log;
            mv %(trim_out)s %(outfile)s;
            ''' % locals()
            outfiles = (outfile,)

        elif self.num_files == 2:
            infile1, infile2 = infiles
            outfile1, outfile2 = outfiles

            cmd = '''trim_galore %(processing_options)s
            --paired
            --phred%(offset)s
            --output_dir %(outdir)s
            %(infile1)s %(infile2)s
            2>>%(output_prefix)s.log;
            mv %(infile1)s_val_1.fq.gz %(outfile1)s;
            mv %(infile2)s_val_2.fq.gz %(outfile2)s;
            ''' % locals()

        return cmd


class Sickle(ProcessTool):

    prefix = "sickle"

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == len(outfiles)
        assert len(infiles) in (1, 2)

        prefix = self.prefix
        offset = Fastq.getOffset("sanger", raises=False)
        processing_options = self.processing_options
        r = {33: 'sanger', 64: 'illumina', 59: 'solexa'}
        quality = r[offset]

        if len(infiles) == 1:
            infile = infiles[0]
            outfile = outfiles[0]
            cmd = '''sickle se
            -g %(processing_options)s
            --qual-type %(quality)s
            --output-file %(outfile)s
            --fastq-file %(infile)s
            2>>%(output_prefix)s.log
            ;''' % locals()

        elif len(infiles) == 2:
            infile1, infile2 = infiles
            outfile1, outfile2 = outfiles
            cmd = '''sickle pe
            -g -s %(processing_options)s
            --qual-type %(quality)s
            -f %(infile1)s -r %(infile2)s
            -o %(outfile1)s -p %(outfile2)s
            2>>%(output_prefix)s.log
            ;''' % locals()

        return cmd


class Trimmomatic(ProcessTool):

    prefix = 'trimmomatic'

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == len(outfiles)
        assert len(infiles) in (1, 2)

        offset = Fastq.getOffset("sanger", raises=False)
        threads = self.threads
        processing_options = self.processing_options
        if len(infiles) == 1:
            infile = infiles[0]
            outfile = outfiles[0]

            cmd = '''trimmomatic SE
            -threads %(threads)i
            -phred%(offset)s
            %(infile)s %(outfile)s
            %(processing_options)s
            2>> %(output_prefix)s.log
            ;''' % locals()

        elif len(infiles) == 2:
            infile1, infile2 = infiles
            outfile1, outfile2 = outfiles

            cmd = '''trimmomatic PE
            -threads %(threads)i
            -phred%(offset)s
            %(infile1)s %(infile2)s
            %(outfile1)s %(output_prefix)s.1.unpaired
            %(outfile2)s %(output_prefix)s.2.unpaired
            %(processing_options)s
            2>> %(output_prefix)s.log;
            checkpoint;
            gzip %(output_prefix)s.*.unpaired;
            ''' % locals()

        return cmd


class FastxTrimmer(ProcessTool):

    prefix = "fastxtrimmer"

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == len(outfiles)
        assert len(infiles) in (1, 2)

        prefix = self.prefix
        offset = Fastq.getOffset("sanger", raises=False)
        processing_options = self.processing_options

        assert len(infiles) == len(outfiles)

        cmds = []
        for infile, outfile in zip(infiles, outfiles):

            cmds.append('''zcat %(infile)s
            | fastx_trimmer
            -Q%(offset)s
            %(processing_options)s
            2>> %(output_prefix)s.log
            | gzip > %(outfile)s
            ;''' % locals())

        return " checkpoint; ".join(cmds)


class Cutadapt(ProcessTool):

    prefix = "cutadapt"

    def build(self, infiles, outfiles, output_prefix):
        prefix = self.prefix
        processing_options = self.processing_options

        assert len(infiles) == len(outfiles)

        cmds = []
        for infile, outfile in zip(infiles, outfiles):

            cmds.append('''zcat %(infile)s
            | cutadapt %(processing_options)s -
            2>> %(output_prefix)s.log
            | gzip > %(outfile)s;''' % locals())

        return " checkpoint; ".join(cmds)


class Reconcile(ProcessTool):

    prefix = "reconcile"

    def process(self, infiles, outfiles, output_prefix):

        assert len(infiles) == 2
        infile1, infile2 = infiles

        cmd = """python %%(scriptsdir)s/fastqs2fastqs.py
        --method=reconcile
        --output-filename-pattern=%(output_prefix)s.fastq.%%s.gz
        %(infile1)s %(infile2)s
        """ % locals()

        return cmd


class Flash(ProcessTool):

    prefix = "flash"

    def get_num_files(self, infiles):
        assert len(infiles) == 2
        return 1

    def process(self, infiles, outfiles, output_prefix):

        prefix = self.prefix
        offset = Fastq.getOffset("sanger", raises=False)
        outdir = os.path.join(os.path.dirname(output_prefix),
                              "flash.dir")
        track = os.path.basename(output_prefix)

        processing_options = self.processing_options

        infile1, infile2 = infiles
        outfile = outfiles[0]

        cmd = '''flash %(infile1)s %(infile2)s
        -p%(offset)s
        %(processing_options)s
        -o %(track)s
        -d %(outdir)s
        2>>%(output_prefix)s.log;
        checkpoint;
        gzip %(outdir)s/*;
        checkpoint;
        mv %(outdir)s/%(track)s.extendedFrags.fastq %(outfile)s;
        ;''' % locals()

        return cmd

    def postprocess(self, infiles):
        # postprocess used to summarise single end fastq twice so that
        # it gets concatenated at the end of the summary tables for
        # both initial paired end files. if a further single end step
        # is included after flash, currently, it will break at the
        # summarise function in pipeline_preprocess

        if self.summarise:
            outdir = self.outdir
            prefix = self.prefix
            infile1 = infiles[0]
            infile_base1 = os.path.basename(infile1)
            infile_base2 = re.sub(".1.fastq.gz", ".2.fastq.gz", infile_base1)
            infile = re.sub(".fastq.1.gz", ".fastq.gz", infile1)
            postprocess_cmd = '''zcat %(infile)s |
            python %%(scriptsdir)s/fastq2summary.py
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base1)s.summary;
            zcat %(infile)s |
            python %%(scriptsdir)s/fastq2summary.py
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base2)s.summary
            ;''' % locals()
        else:
            postprocess_cmd = "checkpoint ;"

        return postprocess_cmd
