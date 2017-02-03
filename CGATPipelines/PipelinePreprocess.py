'''PipelinePreprocess.py - Tasks for processing short read data sets
====================================================================

Processing reads is a common task before mapping. Different sequencing
technologies and applications require different read processing.  This
module provides utility functions to abstract some of these variations

The pipeline does not know what kind of data it gets (a directory
might contain single end or paired end data or both).

The module currently provides modules to perform:
    * hard-trimming (fastx_trimmer, trimmomatic)
    * adapter trimming (trimgalore, trimmomatic, cutadapt)
    * read end quality trimming (trimgalore, trimmomatic)
    * sliding window quality trimming (sickle, trimmomatic)
    * RRBS-specific trimming (trimgalore)

To add a new tool, derive a new class from :class:`ProcessTools`.

Requirements:

* sickle >= 1.33
* fastx >= 0.0.13
* cutadapt >= 1.7.1
* trimmomatic >= 0.32
* trimgalore >= 0.3.3
* flash >= 1.2.6
* pandaseq >= 2.9

Reference
---------

'''

import re
import os
import sqlite3
import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq
import CGATPipelines.PipelineMapping as Mapping
import CGAT.Sra as Sra
import CGAT.Experiment as E


def makeAdaptorFasta(infile, outfile, track, dbh, contaminants_file):
    '''Generate a .fasta file of adaptor sequences that are
    overrepresented in the reads from a sample.

    Requires cutadapt >= 1.7.

    Arguments
    ---------
    infile : string
        Input filename that has been QC'ed. The filename is used to
        check if the input was a :term:`sra` file and guess the
        number of tracks to check.
    outfile : string
        Output filename in :term:`fasta` format.
    track : string
        Track name, used to access FastQC results in database.
    dbh : object
        Database handle.
    contaminants_file : string
        Path of file containing contaminants used for screening by
        Fastqc.

    '''
    tracks = [track]

    if infile.endswith(".sra"):
        # patch for SRA files, look at multiple tracks
        f, fastq_format, datatype = Sra.peek(infile)
        if len(f) == 2:
            tracks = [track + "_fastq_1", track + "_fastq_2"]
    elif infile.endswith(".fastq.1.gz"):
        tracks = [track + "_fastq_1", track + "_fastq_2"]
    elif infile.endswith(".fastq.gz"):
        tracks = [track]

    found_contaminants = []
    for t in tracks:
        table = PipelineTracks.AutoSample(os.path.basename(t)).asTable()

        # if sample name starts with a number, sql table will have
        # prepended "_"
        if re.match("^\d+.*", table):
            table = "_" + table

        query = '''SELECT Possible_Source, Sequence FROM
        %s_fastqc_Overrepresented_sequences;''' % table

        cc = dbh.cursor()
        # if there is no contamination table for even a single sample
        # it will prevent the whole pipeline progressing
        try:
            found_contaminants.extend(cc.execute(query).fetchall())
        except sqlite3.OperationalError:
            E.warn("No table found for {}".format(t))

    if len(found_contaminants) == 0:
        P.touch(outfile)
        return

    # read contaminants from existing file
    with IOTools.openFile(contaminants_file, "r") as inf:
        known_contaminants = [l.split() for l in inf
                              if not l.startswith("#") and l.strip()]
        known_contaminants = {" ".join(x[:-1]): x[-1]
                              for x in known_contaminants}

    # output the full sequence of the contaminant if found
    # in the list of known contaminants, otherwise don't report!

    matched_contaminants = set()
    with IOTools.openFile(outfile, "w") as outf:
        for found_source, found_seq in found_contaminants:
            possible_source = found_source.split(" (")[0]

            if possible_source in known_contaminants:
                matched_contaminants.update((possible_source,))
            else:
                pass

        if len(matched_contaminants) > 0:
            for match in matched_contaminants:
                outf.write(">%s\n%s\n" % (match.replace(" ,", ""),
                                          known_contaminants[match]))


class MasterProcessor(Mapping.SequenceCollectionProcessor):
    """Class for short-read processing tools.

    Processing tools take reads in :term:`fastq` or :term:`sra` formatted
    file and outut one or more :term:`fastq` formatted file.

    This class is a container for processing tools that will
    be applied in sequence. The basic usage is::

        m = MasterProcessor()
        m.add(FastxTrimmer())
        m.add(Trimmomatic())

        statement = m.build(("Sample1.fastq.gz", ),
                            "processed.dir/",
                            "Sample1")
        Pipeline.run()

    Attributes
    ----------
    save : bool
        If True, save intermediate files in the final output
        directory. Otherwise, a temporary directory is used.
    summarize : bool
        If True, summarize the output after each step.
    threads : int
        Number of processing threads to use.
    outdir : string
        Output directory for intermediate files.
    processors : list
        List of processors that will be applied in sequence.

    Arguments
    ---------
    save : bool
        If True, save intermediate files in the final output
        directory. Otherwise, a temporary directory is used.
    summarize : bool
        If True, summarize the output after each step.
    threads : int
        Number of processing threads to use.

    """

    # compress temporary fastq files with gzip
    compress = False

    def __init__(self,
                 save=True,
                 summarize=False,
                 threads=1,
                 qual_format='phred64',
                 *args, **kwargs):
        self.save = save
        self.summarize = summarize
        self.threads = threads
        if self.save:
            self.outdir = "processed.dir"
        else:
            self.outdir = P.getTempDir(shared=True)

        self.processors = []
        self.qual_format = qual_format

    def add(self, processor):
        """add a processor to the list of tools to be executed."""
        self.processors.append(processor)

    def cleanup(self):
        '''return statement to clean up temporary files after processing.'''
        if self.save:
            statement = 'checkpoint;'
        else:
            statement = 'rm -rf %s;' % self.outdir
        statement += '''rm -rf %s;''' % (self.tmpdir_fastq)
        return statement

    def build(self, infile, output_prefix, track):
        '''build a preprocessing statement.

        This method iterates through all processing tools
        that have been added and builds a single command line
        statement executing them in sequence.

        Arguments
        ---------
        infile : string
            Input filename.
        output_prefix : string
            Output prefix to use for files created by a processing
            statement.
        track : string
            Prefix for output files. The suffix ".fastq.gz" will
            be appended for single-end data sets. For paired end
            data sets, two files will be created ending in
            ".fastq.1.gz" and ".fastq.2.gz"

        Returns
        -------
        statement : string
            A command line statement to be executed in a pipeline.

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

            if idx == len(self.processors) - 1:
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
                        | cgat fastq2summary
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
    '''defines class attributes for a sequence utility tool

    Derived classes need to implement a :func:`ProcessTools.build`
    method.

    Attributes
    ----------
    processing_options : string
        List of options to send to the processing tool. These
        are tool specific.
    threads : int
        Number of threads to use.

    Arguments
    ---------
    options : string
        List of options to send to the processing tool. These
        are tool specific.
    threads : int
        Number of threads to use.
    '''

    def __init__(self, options, threads=1, untrimmed=0):
        self.processing_options = options
        self.threads = threads
        self.untrimmed = untrimmed

    def get_num_files(self, infiles):
        """return the number of outputfiles created by
        processing `infiles`"""

        # the default is that the same number of files are
        # returned
        return len(infiles)

    def build(self, infiles, outfiles, output_prefix):
        """build a command line statement executing the tool.

        Arguments
        ---------
        infiles : list
           List of input filenames in :term:`fastq` format.
        outfiles : list
           List of output filenames in :term:`fastq` format.
        output_prefix : list
           Prefix to add to PATH of additional output files created by
           the tool.

        .. note::
            This method needs to overloaded by derived classes.

        """

        raise NotImplementedError(
            "build() method needs to be implemented in derived class")


class Trimgalore(ProcessTool):
    """Read processing - run trimgalore"""

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
    """Read processing - run sickle"""

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
    """Read processing - run trimmomatic"""

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
    """Read processing - run fastxtrimmer"""

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
    """Read processing - run cutadapt"""

    prefix = "cutadapt"

    def __init__(self, options, threads=1, process_paired=0, *args, **kwargs):
        self.process_paired = process_paired
        ProcessTool.__init__(self, options, threads, *args, **kwargs)

    def build(self, infiles, outfiles, output_prefix):
        prefix = self.prefix
        processing_options = self.processing_options
        untrimmed = self.untrimmed

        assert len(infiles) == len(outfiles)

        cmds = []
        if self.process_paired and len(infiles) == 2:
            in1, in2 = infiles
            out1, out2 = outfiles

            if "fastq" in in1:
                format = "--format=fastq"
            elif "fasta" in in1:
                format = "--format=fasta"
            else:
                format = ""

            untrimmed_output1, untrimmed_output2 = \
                [i.replace(".fast", "_untrimmed.fast")
                 for i in infiles]

            if untrimmed:
                processing_options += \
                    "--untrimmed-output=%(untrimmed_output1)s" \
                    "--untrimmed-output-paired=%(untrimmed_output2)s" % locals()

            cmds.append('''
            cutadapt %(processing_options)s %(in1)s %(in2)s
                     -p %(out2)s -o %(out1)s %(format)s
            2>> %(output_prefix)s.log; ''' % locals())

            if untrimmed:
                cmds.append("gzip %s;" % untrimmed_output1)
                cmds.append("gzip %s;" % untrimmed_output2)

        else:
            for infile, outfile in zip(infiles, outfiles):
                outfile_untrimmed = outfile.replace(".fastq",
                                                    "_untrimmed.fastq")
                if untrimmed:
                    processing_options += " --untrimmed-output=%s" % \
                        outfile_untrimmed

                cmds.append('''zcat %(infile)s
                | cutadapt %(processing_options)s -
                2>> %(output_prefix)s.log
                | gzip > %(outfile)s;''' % locals())

                if untrimmed:
                    cmds.append("gzip %s;" % outfile_untrimmed)

        return " checkpoint; ".join(cmds)


class Reconcile(ProcessTool):
    """Read processing - reconcile read names in two :term:`fastq` files.

    This method iterates through two :term:`fastq` files and removes
    those entries that do not appear in both.

    The :term:`fastq` files are expected to be sorted in the same
    order.

    """

    prefix = "reconcile"

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == 2
        infile1, infile2 = infiles

        cmd = """cgat fastqs2fastqs
        --method=reconcile
        --output-filename-pattern=%(output_prefix)s.fastq.%%%%s.gz
        %(infile1)s %(infile2)s;
        """ % locals()

        return cmd


class Flash(ProcessTool):
    """Read processing - run flash"""

    prefix = "flash"

    def get_num_files(self, infiles):
        assert len(infiles) == 2
        return 1

    def build(self, infiles, outfiles, output_prefix):

        prefix = self.prefix
        offset = Fastq.getOffset("sanger", raises=False)
        outdir = os.path.join(output_prefix + ".dir")
        track = os.path.basename(output_prefix)

        processing_options = self.processing_options

        infile1, infile2 = infiles
        outfile = outfiles[0]

        cmd = '''flash %(infile1)s %(infile2)s
        -p %(offset)s
        %(processing_options)s
        -o %(track)s
        -d %(outdir)s
        >& %(output_prefix)s-flash.log;
        checkpoint;
        gzip %(outdir)s/*;
        checkpoint;
        mv %(outdir)s/%(track)s.extendedFrags.fastq.gz %(outfile)s;
        ''' % locals()

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
            cgat fastq2summary
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base1)s.summary;
            zcat %(infile)s |
            cgat fastq2summary
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base2)s.summary
            ;''' % locals()
        else:
            postprocess_cmd = "checkpoint ;"

        return postprocess_cmd


class ReverseComplement(ProcessTool):
    """Reverse complement fastq files."""

    prefix = "reversecomplement"

    def build(self, infiles, outfiles, output_prefix):

        assert len(infiles) == len(outfiles)

        cmds = []
        for infile, outfile in zip(infiles, outfiles):
            cmds.append('''zcat %(infile)s
            | cgat fastq2fastq
            --method=reverse-complement
            --log=%(output_prefix)s.log
            | gzip > %(outfile)s;
            ''' % locals())

        return " checkpoint; ".join(cmds)


class Pandaseq(ProcessTool):
    """Read processing - run pandaseq"""

    prefix = "pandaseq"

    def get_num_files(self, infiles):
        print("infiles=", infiles)
        assert len(infiles) == 2
        return 1

    def build(self, infiles, outfiles, output_prefix):

        prefix = self.prefix
        offset = Fastq.getOffset("sanger", raises=False)
        outdir = os.path.join(output_prefix + ".dir")
        track = os.path.basename(output_prefix)

        processing_options = self.processing_options
        threads = self.threads

        infile1, infile2 = infiles
        outfile = outfiles[0]

        cmd = '''pandaseq -f %(infile1)s -r %(infile2)s
        %(processing_options)s
        -T %(threads)i
        -U >(gzip > %(outfile)s.unpaired.gz)
        -w >(gzip > %(outfile)s)
        -F
        -G %(output_prefix)s-pandaseq.log.bgz;
        >& %(output_prefix)s-pandaseq.log;
        checkpoint;
        gzip %(outdir)s/*;
        checkpoint;
        ''' % locals()

        return cmd
