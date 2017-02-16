'''PipelineMapping.py - Utility functions for mapping short reads
==============================================================

:Release: $Id$
:Date: |today|
:Tags: Python

Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra`
files) of different data (single end, paired end) with different
mapping algorithms (bowtie, tophat, stampy). This module provides
utility functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra`
archive might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

It implements:

   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

The basic class :class:`SequenceCollectionProcessor`.  The aim of this
class and its derivatives is to build a sequence of command line
statements that can be send to a single node on a cluster to process
the input data.

The basic usage inside a pipeline task is as such::

    @transform()
    def mapReads(infile, outfile):

        # initialize the Tool
        m = PipelineMapping.Hisat(
             executable=P.substituteParameters(**locals())["hisat_executable"],
             strip_sequence=PARAMS["strip_sequence"])

        # build the command line statement
        statement = m.build((infile,), outfile)

        P.run()

When implementing a tool, avoid specifying algorithmic options as
class variables. Instead use an option string that can be set in
:file:`pipeline.ini`. The only arguments to a tool constructor should
pertain to pipeline integration, such as filenames, index locations,
threading and in general processing options that change the tools
input/output, as these need to be tracked by the pipeline.

The benefit of this approach is to provide compleeat control to the
user and is likely to work with different versions of a tool, i.e., if
a command line option to a tool changes, just the configuration file
needs to be changed, but the code remains the same.

Requirements:

* cufflinks >= 2.2.1
* fastq-dump >= 2.1.7
* fastqc >= 0.11.2
* fastq_screen >= 0.4.4
* sailfish >= 0.6.3
* picardtools >= 1.106
* samtools >= 1.1
* tophat >= 2.0.13 (optional)
* bowtie >= 1.0.0 (optional)
* bowtie2 >= 2.2.3 (optional)
* bwa >= 0.7.8
* gsnap >= 2014-01-21 (optional)
* star >= 2.3.0e (optional)
* bismark >= 0.12.5 (optional)
* stampy >= 1.0.23 (optional)
* butter >= 0.3.2 (optional)
* hisat >= 0.1.4 (optional)
* shortstack >3.4 (optional)

Reference
---------

'''

import os
import glob
import collections
import re
import gzip
import itertools
import CGATPipelines.Pipeline as P
import logging as L
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Fastq as Fastq
import CGAT.IndexedFasta as IndexedFasta
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGAT.Sra as Sra

SequenceInformation = collections.namedtuple("SequenceInformation",
                                             """paired_end
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""")


def getSequencingInformation(track):
    '''glean sequencing information from *track*.

    This method looks for possible filename suffixes for
    track that exist. Depending on the files present it
    will return a named tuple with the following columns:

    paired_end
       True if data is paired end
    filename_first
       Filename of first read in pair
    filename_second
       Filename of second read in pair (if paired)
    readlength_first
       Read length of first read in pair
    readlength_second
       Read length of second read in pair (if given)

    Arguments
    ---------
    track : string
        A track name

    Returns
    -------
    data : tuple
    '''

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
        second_length = Fastq.getReadLength(second_pair)

    return SequenceInformation._make((second_pair is not None,
                                      first_pair, second_pair,
                                      Fastq.getReadLength(first_pair),
                                      second_length,
                                      colour))


def mergeAndFilterGTF(infile, outfile, logfile,
                      genome,
                      max_intron_size=None,
                      remove_contigs=None,
                      rna_file=None):
    '''sanitize transcripts file for cufflinks analysis.

    Merge exons separated by small introns (< 5bp).

    Transcripts will be ignored that
       * have very long introns (`max_intron_size`) (otherwise,
         cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)

    Optionally remove transcripts based on repetitive sequences by
       supplying a repetitve `rna_file`.

    This method preserves all features in a gtf file (exon, CDS, ...).

    Arguments
    ---------
    infile : string
       Input filename in :term:`gtf` format
#    outfile : string
       Output filename in :term:`gtf` format
    logfile : string
       Output filename for logging information.
    genome : string
       Filename (without extension) of indexed genome file
       in :term:`fasta` format.
    max_intron_size : int
       Remove transripts with introns larger than this value.
    remove_contigs : string
       Remove transcripts on contigs matching this regular
       expression string.
    rna_file : string
       Filename of :term:`gff` formatted file with repetetive
       sequences. If given, all transcripts overlapping any regions
       in this file will be removed.

    Returns
    -------
    kept_genes : dict
        a dictionary of all gene_ids that have been kept.

    '''

    c = E.Counter()

    outf = gzip.open(outfile, "w")

    E.info("filtering by contig and removing long introns")
    contigs = set(IndexedFasta.IndexedFasta(genome).getContigs())

    rx_contigs = None
    #
    if remove_contigs is not None:
        rx_contigs = re.compile(remove_contigs)
        E.info("removing contigs %s" % remove_contigs)

    rna_index = None
    if rna_file is not None:
        if not os.path.exists(rna_file):
            E.warn("file '%s' to remove repetetive rna does not exist" %
                   rna_file)
        else:
            rna_index = GTF.readAndIndex(
                GTF.iterator(IOTools.openFile(rna_file, "r")))
            E.info("removing ribosomal RNA in %s" % rna_file)

    gene_ids = {}

    logf = IOTools.openFile(logfile, "w")
    logf.write("gene_id\ttranscript_id\treason\n")

    for all_exons in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile))):

        c.input += 1

        e = all_exons[0]
        # filtering
        if e.contig not in contigs:
            c.missing_contig += 1
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "missing_contig")) + "\n")
            continue

        if rx_contigs and rx_contigs.search(e.contig):
            c.remove_contig += 1
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "remove_contig")) + "\n")
            continue

        if rna_index and all_exons[0].source != 'protein_coding':
            found = False
            for exon in all_exons:
                if rna_index.contains(e.contig, e.start, e.end):
                    found = True
                    break
            if found:
                logf.write(
                    "\t".join((e.gene_id, e.transcript_id,
                               "overlap_rna")) + "\n")
                c.overlap_rna += 1
                continue

        is_ok = True

        # keep exons and cds separate by grouping by feature
        all_exons.sort(key=lambda x: x.feature)
        new_exons = []

        for feature, exons in itertools.groupby(
                all_exons, lambda x: x.feature):

            tmp = sorted(list(exons), key=lambda x: x.start)

            gene_ids[tmp[0].transcript_id] = tmp[0].gene_id

            l, n = tmp[0], []

            for e in tmp[1:]:
                d = e.start - l.end
                if max_intron_size and d > max_intron_size:
                    is_ok = False
                    break
                elif d < 5:
                    l.end = max(e.end, l.end)
                    c.merged += 1
                    continue

                n.append(l)
                l = e

            n.append(l)
            new_exons.extend(n)

            if not is_ok:
                break

        if not is_ok:
            logf.write(
                "\t".join((e.gene_id, e.transcript_id,
                           "bad_transcript")) + "\n")
            c.skipped += 1
            continue

        new_exons.sort(key=lambda x: (x.start, x.gene_id, x.transcript_id))

        for e in new_exons:
            outf.write("%s\n" % str(e))
            c.exons += 1

        c.output += 1

    outf.close()
    L.info("%s" % str(c))

    return gene_ids


def resetGTFAttributes(infile, genome, gene_ids, outfile):
    """set GTF attributes in :term:`gtf` formatted file so that they are
    compatible with cufflinks.

    This method runs cuffcompare with `infile` against itself to add
    attributes such as p_id and tss_id.

    Arguments
    ---------
    infile : string
        Filename of :term:`gtf`-formatted input file
    genome : string
       Filename (without extension) of indexed genome file
       in :term:`fasta` format.
    gene_ids : dict
       Dictionary mapping transcript ids to gene ids.
    outfile : string
       Output filename in :term:`gtf` format
    """
    tmpfile1 = P.getTempFilename(".")
    tmpfile2 = P.getTempFilename(".")

    #################################################
    E.info("adding tss_id and p_id")

    # The p_id attribute is set if the fasta sequence is given.
    # However, there might be some errors in cuffdiff downstream:
    #
    # cuffdiff: bundles.cpp:479: static void HitBundle::combine(const std::
    # vector<HitBundle*, std::allocator<HitBundle*> >&, HitBundle&): Assertion
    # `in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id()' failed.
    #
    # I was not able to resolve this, it was a complex
    # bug dependent on both the read libraries and the input reference gtf
    # files
    job_memory = "5G"

    statement = '''
    cuffcompare -r <( gunzip < %(infile)s )
         -T
         -s %(genome)s.fa
         -o %(tmpfile1)s
         <( gunzip < %(infile)s )
         <( gunzip < %(infile)s )
    > %(outfile)s.log
    '''
    P.run()

    #################################################
    E.info("resetting gene_id and transcript_id")

    # reset gene_id and transcript_id to ENSEMBL ids
    # cufflinks patch:
    # make tss_id and p_id unique for each gene id
    outf = IOTools.openFile(tmpfile2, "w")
    map_tss2gene, map_pid2gene = {}, {}
    inf = IOTools.openFile(tmpfile1 + ".combined.gtf")

    def _map(gtf, key, val, m):
        if val in m:
            while gene_id != m[val]:
                val += "a"
                if val not in m:
                    break
        m[val] = gene_id

        gtf.setAttribute(key, val)

    for gtf in GTF.iterator(inf):
        transcript_id = gtf.oId
        gene_id = gene_ids[transcript_id]
        gtf.setAttribute("transcript_id", transcript_id)
        gtf.setAttribute("gene_id", gene_id)

        # set tss_id
        try:
            tss_id = gtf.tss_id
        except AttributeError:
            tss_id = None
        try:
            p_id = gtf.p_id
        except AttributeError:
            p_id = None

        if tss_id:
            _map(gtf, "tss_id", tss_id, map_tss2gene)
        if p_id:
            _map(gtf, "p_id", p_id, map_pid2gene)

        outf.write(str(gtf) + "\n")

    outf.close()

    # sort gtf file
    PipelineGeneset.sortGTF(tmpfile2, outfile)

    # make sure tmpfile1 is NEVER empty
    assert tmpfile1
    for x in glob.glob(tmpfile1 + "*"):
        os.unlink(x)
    os.unlink(tmpfile2)


class SequenceCollectionProcessor(object):
    """base class for processors of sequence collections.

    Processors of sequence collections are mappers, trimmers, etc.

    The aim of these processors is to build a sequence
    of command line statements that can be send to a single node
    on a cluster to process the input data.

    Attributes
    ==========
    compress : bool
        If True, compress temporary fastq files with gzip
    convert : bool
        If True, convert quality scores to Sanger quality scores.
    preserve_colourspace : bool
        If True, preserve colour space files. By default they are
        converted to fastq.
    tmpdir_fastq : string
        Directory with the locations of temporary :term:`fastq`
        formatted files. This directory can be used as a general
        temporary directory by a mapper.
    """

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    # set to True if you want to preserve colour space files.
    # By default, they are converted to fastq.
    preserve_colourspace = False

    # Temporary directory where temporary fastq files will
    # be located.
    tmpdir_fastq = None

    keep_sra = True

    def __init__(self, *args, **kwargs):
        pass

    def quoteFile(self, filename):
        '''return a quoted file for in-situ uncompression.

        The quoting adds uncompression for compressed files for
        programs that expect uncompressed files.

        .. note::
            This will only work if the downstream programs reads
            the file only once.

        '''
        if filename.endswith(".gz") and not self.compress:
            return "<( gunzip < %s )" % filename
        else:
            return filename

    def preprocess(self, infiles, outfile):
        '''build a preprocessing statement

        Build a command line statement that extracts/converts various
        input formats to :term:`fastq` formatted files.

        Mapping qualities are changed to solexa format.

        Returns
        -------
        statement : string
            The command line statement for pre-processing.
        files : list
            List of :term:`fastq` formatted files that will
            be created by `statement`.
        '''
        # explicitly assign qual_format to use in string formatting
        try:
            assert self.qual_format
            qual_format = self.qual_format
        except AttributeError:
            qual_format = 'illumina-1.8'

        assert len(infiles) > 0, "no input files for processing"

        tmpdir_fastq = P.getTempDir(shared=True)

        # create temporary directory again for nodes
        statement = ["mkdir -p %s" % tmpdir_fastq]
        fastqfiles = []

        # get track by extension of outfile
        track = os.path.splitext(os.path.basename(outfile))[0]

        if self.compress:
            compress_cmd = "| gzip"
            extension = ".gz"
        else:
            compress_cmd = ""
            extension = ""

        for infile in infiles:

            if infile.endswith(".export.txt.gz"):
                # single end illumina export
                statement.append("""gunzip < %(infile)s
                | awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
                {if ($1 != "")
                {readname=sprintf("%%%%s_%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $2, $3, $4, $5, $6);}
                else {readname=sprintf("%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $3, $4, $5, $6);}
                printf("@%%%%s\\n%%%%s\\n+\\n%%%%s\\n",readname,$9,$10);}'
                %(compress_cmd)s
                > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals())
                fastqfiles.append(
                    ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))
            elif infile.endswith(".fa.gz"):
                statement.append(
                    '''gunzip < %(infile)s %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fa''' % locals())
                fastqfiles.append(("%s/%s.fa" % (tmpdir_fastq, track),))
                self.datatype = "fasta"

            elif infile.endswith(".remote"):
                files = []
                for line in IOTools.openFile(infile):
                    repo, acc = line.strip().split("\t")[:2]
                    if repo == "SRA":
                        statement.append(Sra.prefetch(acc))
                        f, format = Sra.peek(acc)
                        statement.append(Sra.extract(acc, tmpdir_fastq))

                        extracted_files = ["%s/%s" % (
                            tmpdir_fastq, os.path.basename(x))
                            for x in sorted(f)]
                        if not self.keep_sra:
                            statement.append(Sra.clean_cache(acc))
                        files.extend(extracted_files)

                    elif repo == "ENA":
                        filenames, dl_paths = Sra.fetch_ENA_files(acc)
                        for f in dl_paths:
                            statement.append(Sra.fetch_ENA(f, tmpdir_fastq))
                        files.extend([os.path.join(tmpdir_fastq, x) for x
                                      in filenames])

                    elif repo == "TCGA":
                        tar_name = line.strip().split("\t")[2]
                        token = glob.glob("gdc-user-token*")

                        if len(token) > 0:
                            token = token[0]
                        else:
                            token = None

                        statement.append(Sra.fetch_TCGA_fastq(acc,
                                                              tar_name,
                                                              token,
                                                              tmpdir_fastq))

                        files.append(os.path.join(
                            tmpdir_fastq, acc + "_1.fastq.gz"))
                        files.append(os.path.join(
                            tmpdir_fastq, acc + "_2.fastq.gz"))

                    else:
                        raise ValueError("Unknown repository: %s" % repo)

                if len(files) == 1 and files[0].endswith("_1.fastq.gz"):
                    new_files = [re.sub("_1.fastq.gz",
                                        ".fastq.gz",
                                        files[0])]

                elif len(files) == 2:
                    new_files = [re.sub(r"_([12]).fastq.gz",
                                        r".fastq.\1.gz",
                                        x) for x in files]

                else:
                    new_files = files

                out_base = os.path.basename(os.path.splitext(outfile)[0])
                new_files = [re.sub(r"%s/(.+).(fastq..*gz)" % tmpdir_fastq,
                                    r"%s/%s_\1.\2" % (tmpdir_fastq, out_base),
                                    nf) for nf in new_files]

                for old, new in zip(files, new_files):
                    if old != new:
                        statement.append("mv %s %s" % (old, new))

                fastqfiles.append(new_files)

            elif infile.endswith(".sra"):
                # sneak preview to determine if paired end or single end
                outdir = P.getTempDir()
                f, format, datatype = Sra.peek(infile)
                E.info("sra file contains the following files: %s" % f)

                # T.S need to use abi-dump for colorspace files
                if datatype == "basecalls":
                    tool = "fastq-dump"
                    self.datatype = "basecalls"
                elif datatype == "colorspace":
                    tool = "abi-dump"
                    self.datatype = "solid"

                # add extraction command to statement
                statement.append(Sra.extract(infile, tmpdir_fastq, tool))

                sra_extraction_files = ["%s/%s" % (
                    tmpdir_fastq, os.path.basename(x)) for x in sorted(f)]

                # if format is not sanger, add convert command to statement
                if 'sanger' in format and self.convert:

                    # single end fastq
                    if len(sra_extraction_files) == 1:

                        infile = sra_extraction_files[0]
                        track = os.path.splitext(os.path.basename(infile))[0]
                        if track.endswith("_1.fastq"):
                            track = track[:-8]
                        statement.append("""gunzip < %(infile)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log
                        %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.fastq%(extension)s
                        """ % locals())

                        fastqfiles.append(
                            ("%(tmpdir_fastq)s/%(track)s_converted.fastq%(extension)s"
                             % locals(),))

                    # paired end fastqs
                    elif len(sra_extraction_files) == 2:

                        infile, infile2 = sra_extraction_files
                        track = os.path.splitext(os.path.basename(infile))[0]

                        statement.append("""gunzip < %(infile)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.1.fastq%(extension)s;
                        gunzip < %(infile2)s
                        | cgat fastq2fastq
                        --method=change-format --target-format=sanger
                        --guess-format=phred64
                        --log=%(outfile)s.log %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s_converted.2.fastq%(extension)s
                        """ % locals())

                        fastqfiles.append(
                            ("%s/%s_converted.1.fastq%s" %
                             (tmpdir_fastq, track, extension),
                             "%s/%s_converted.2.fastq%s" %
                             (tmpdir_fastq, track, extension)))
                else:
                    # Usually sra extraction files have format
                    # '1_fastq.gz' for both single and paired end files
                    # This code corrects the output to the format expected by
                    # CGAT Pipelines
                    infile = sra_extraction_files[0]
                    basename = os.path.basename(infile)

                    if(len(sra_extraction_files) == 1 and
                       basename.endswith("_1.fastq.gz") and
                       self.datatype != "solid"):

                        basename = basename[:-11] + ".fastq.gz"
                        statement.append(
                            "mv %s %s/%s" % (infile, tmpdir_fastq, basename))
                        fastqfiles.append(
                            ("%s/%s" % (tmpdir_fastq, basename),))

                    # if sra extracted file is SOLID, need to keep
                    # record of qual files
                    elif self.datatype == "solid":
                        # single end SOLiD data
                        infile = P.snip(infile, "_1.fastq.gz") + \
                            "_F3.csfasta.gz"
                        quality = P.snip(
                            infile, "_F3.csfasta.gz") + "_F3_QV.qual.gz"

                        # qual file does not exist as tmpdir from
                        # SRA.extract is removed
                        # if not os.path.exists(quality):
                        #    raise ValueError("no quality file for %s" % infile)

                        fastqfiles.append((infile, quality))

                    # T.S I'm not sure if this works. Need a test case!
                    elif infile.endswith(".csfasta.F3.gz"):
                        # paired end SOLiD data
                        bn = P.snip(infile, ".csfasta.F3.gz")
                        # order is important - mirrors tophat reads followed by
                        # quals
                        f = []
                        for suffix in ("csfasta.F3", "csfasta.F5",
                                       "qual.F3", "qual.F5"):
                            fn = "%(bn)s.%(suffix)s" % locals()
                            if not os.path.exists(fn + ".gz"):
                                raise ValueError(
                                    "expected file %s.gz missing" % fn)
                            statement.append("""gunzip < %(fn)s.gz
                            %(compress_cmd)s
                            > %(tmpdir_fastq)s/%(basename)s.%(suffix)s%(extension)s
                            """ %
                                             locals())
                            f.append(
                                "%(tmpdir_fastq)s/%(basename)s.%(suffix)s%(extension)s" %
                                locals())
                        fastqfiles.append(f)

                    elif len(sra_extraction_files) == 2:
                        infile2 = sra_extraction_files[1]
                        if basename.endswith("_1.fastq.gz"):
                            basename1 = basename[:-11] + ".fastq.1.gz"
                            basename2 = basename[:-11] + ".fastq.2.gz"
                        statement.append(
                            "mv %s %s/%s; mv %s %s/%s" %
                            (infile, tmpdir_fastq, basename1,
                             infile2, tmpdir_fastq, basename2))
                        fastqfiles.append(
                            ("%s/%s" % (tmpdir_fastq, basename1),
                             "%s/%s" % (tmpdir_fastq, basename2)))

                    else:
                        fastqfiles.append(sra_extraction_files)

            elif infile.endswith(".fastq.gz"):
                format = Fastq.guessFormat(
                    IOTools.openFile(infile, "r"), raises=False)
                if 'sanger' not in format and self.convert:
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))
                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.append((infile, ))

            elif infile.endswith(".csfasta.gz"):
                # single end SOLiD data
                if self.preserve_colourspace:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"
                    if not os.path.exists(quality):
                        raise ValueError("no quality file for %s" % infile)
                    statement.append("""gunzip < %(infile)s
                    > %(tmpdir_fastq)s/%(track)s.csfasta%(extension)s""" %
                                     locals())
                    statement.append("""gunzip < %(quality)s
                    > %(tmpdir_fastq)s/%(track)s.qual%(extension)s""" %
                                     locals())
                    fastqfiles.append(("%s/%s.csfasta%s" %
                                       (tmpdir_fastq, track, extension),
                                       "%s/%s.qual%s" %
                                       (tmpdir_fastq, track, extension)))
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))

            elif infile.endswith(".csfasta.F3.gz"):
                # paired end SOLiD data
                if self.preserve_colourspace:
                    bn = P.snip(infile, ".csfasta.F3.gz")
                    # order is important - mirrors tophat reads followed by
                    # quals
                    f = []
                    for suffix in ("csfasta.F3", "csfasta.F5",
                                   "qual.F3", "qual.F5"):
                        fn = "%(bn)s.%(suffix)s" % locals()
                        if not os.path.exists(fn + ".gz"):
                            raise ValueError(
                                "expected file %s.gz missing" % fn)
                        statement.append("""gunzip < %(fn)s.gz
                        %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s
                        """ %
                                         locals())
                        f.append(
                            "%(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s" %
                            locals())
                    fastqfiles.append(f)
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))

            elif infile.endswith(".fastq.1.gz"):

                bn = P.snip(infile, ".fastq.1.gz")
                infile2 = "%s.fastq.2.gz" % bn
                if not os.path.exists(infile2):
                    raise ValueError(
                        "can not find paired ended file "
                        "'%s' for '%s'" % (infile2, infile))

                format = Fastq.guessFormat(
                    IOTools.openFile(infile), raises=False)

                if 'sanger' not in format and qual_format != 'phred64':
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                    checkpoint;
                    gunzip < %(infile2)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                    """ % locals())
                    fastqfiles.append(
                        ("%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                         "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension)))

                elif 'sanger' not in format and qual_format == 'phred64':
                    statement.append("""gunzip < %(infile)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                    checkpoint;
                    gunzip < %(infile2)s
                    | cgat fastq2fastq
                    --method=change-format --target-format=sanger
                    --guess-format=%(qual_format)s
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                    """ % locals())
                    fastqfiles.append(
                        ("%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                         "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension)))
                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.append((infile,
                                       infile2,))

            else:
                raise NotImplementedError("unknown file format %s" % infile)

        self.tmpdir_fastq = tmpdir_fastq

        assert len(fastqfiles) > 0, "no fastq files for mapping"
        return "; checkpoint;".join(statement) + ";", fastqfiles


class Mapper(SequenceCollectionProcessor):
    '''Base class for short-read mappers.

    These tools map reads in :term:`fastq` or :term:`sra` formatted
    files and outut a :term:`BAM` formatted file.

    .. note::
        This class just sets the attributes below, but is not
        implementing them.

    Attributes
    ----------
    datatype : string
        Datatype of input
    strip_sequence : bool
        If True, remove sequence from BAM files in post-processing
        step.
    remove_non_unique : bool
        If True, remove non-unique matches from the BAM file. This
        processing will happen in addition to any filters applied
        in the mapper.

    Arguments
    ---------
    executable : string
        Executable to use. If unset, use default.
    strip_sequence : bool
        If True, remove sequence from BAM files in post-processing.
    remove_non_unique : bool
        If True, remove non-unique matches from the BAM file. This
        processing will happen in addition to any filters applied
        in the mapper.
    tool_options : string
        Options to be passed to the processing tool.

    '''

    datatype = "fastq"

    # strip bam files of sequenca and quality information
    strip_sequence = False

    # remove non-unique matches in a post-processing step.
    # Many aligners offer this option in the mapping stage
    # If only unique matches are required, it is better to
    # configure the aligner as removing in post-processing
    # adds to processing time.
    remove_non_unique = False

    def __init__(self,
                 executable=None,
                 strip_sequence=False,
                 remove_non_unique=False,
                 tool_options="",
                 *args, **kwargs):
        '''
       '''
        SequenceCollectionProcessor.__init__(self, *args, **kwargs)

        if executable:
            self.executable = executable
        self.strip_sequence = strip_sequence
        self.remove_non_unique = remove_non_unique

        # tool options to be passed on to the mapping tool
        self.tool_options = tool_options

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles.
        '''
        return ""

    def postprocess(self, infiles, outfile):
        '''collect output data and postprocess.'''
        return ""

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s;''' % (self.tmpdir_fastq)

        return statement

    def build(self, infiles, outfile):
        '''run mapper

        This method combines the output of the :meth:`preprocess`,
        :meth:`mapper`, :meth:`postprocess` and :meth:`clean` sections
        into a single statement.

        Arguments
        ---------
        infiles : list
             List of input filenames
        outfile : string
             Output filename

        Returns
        -------
        statement : string
             A command line statement. The statement can be a series
             of commands separated by ``;`` and/or can be unix pipes.

        '''

        cmd_preprocess, mapfiles = self.preprocess(infiles, outfile)
        cmd_mapper = self.mapper(mapfiles, outfile)
        cmd_postprocess = self.postprocess(infiles, outfile)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";"),\
            "missing ';' at end of command %s" % cmd_preprocess.strip()
        assert cmd_mapper.strip().endswith(";"),\
            "missing ';' at end of command %s" % cmd_mapper.strip()
        if cmd_postprocess:
            assert cmd_postprocess.strip().endswith(";"),\
                "missing ';' at end of command %s" % cmd_postprocess.strip()
        if cmd_clean:
            assert cmd_clean.strip().endswith(";"),\
                "missing ';' at end of command %s" % cmd_clean.strip()

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_mapper,
                                          cmd_postprocess,
                                          cmd_clean))

        return statement


class FastQc(Mapper):
    """run the FastQC_ tool.

    Arguments
    ---------
    outdir : string
        Output directory
    contaminants : string
        Filename with sequence for contaminant sequences to check for.

    """

    compress = True

    def __init__(self,
                 outdir=".",
                 contaminants=None,
                 qual_format='phred64',
                 *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.outdir = outdir
        if contaminants:
            self.contaminants = self.parse_contaminants(contaminants)
        else:
            self.contaminants = None
        self.qual_format = qual_format

    def parse_contaminants(self, contaminants):
        '''remove misplaced tabs from contaminants file.

        Misplated tab characters in the contaminants file can cause
        Fastqc to hiccup when provided as an adaptor contaminants
        file.  Parse file and insert a single tab.

        Using this file is only appropriate if the adaptors sequences
        are shorter than the reads.
        '''

        # if no contaminants file provided return None
        try:
            assert contaminants
        except:
            return None

        try:
            assert os.path.exists(contaminants)
        except AssertionError:
            raise Exception(contaminants)

        # read in file and split into adaptor/sequence
        adaptor_dict = {}
        with IOTools.openFile(contaminants, "r") as ofile:
            for line in ofile.readlines():
                if line[0] == '#':
                    pass
                else:
                    adapt = line.split("\t")
                    adaptor = adapt[0]
                    sequence = adapt[-1]
                    adaptor_dict[adaptor] = sequence
        ofile.close()

        # get temporary file name
        outfile = P.getTempFilename(shared=True)
        with IOTools.openFile(outfile, "w") as wfile:
            for key, value in list(adaptor_dict.items()):
                wfile.write("%s\t%s\n" % (key, value))
        wfile.close()

        return outfile

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles.

        The output is created in outdir.

        The FASTQC output files are extracted in order to available in
        the pipeline report.

        '''

        contaminants = self.contaminants
        outdir = self.outdir
        statement = []
        for f in infiles:
            for i, x in enumerate(f):
                track = os.path.basename(re.sub(".fastq.*", "", x))

                if contaminants:
                    contaminants_cmd = "-a %s" % contaminants
                else:
                    contaminants_cmd = ""

                statement.append(
                    '''fastqc --extract --outdir=%(outdir)s %(x)s
                    %(contaminants_cmd)s >& %(outfile)s ; ''' % locals())
        statement.append('''rm -f %(contaminants)s ;''' % locals())
        return " ".join(statement)


class FastqScreen(Mapper):
    '''run Fastq_screen to test contamination by other organisms.'''

    compress = True

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles
        The output is created in outdir. The output files
        are extracted.
        '''

        if len(infiles) > 1:
            raise ValueError(
                "Fastq_screen can only operate on one fastq file at a time")

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        if nfiles == 1:
            input_files = infiles[0][0]
        elif nfiles == 2:
            infile1, infile2 = infiles[0]
            input_files = '''-- paired %(infile1)s %(infile2)s''' % locals()
        else:
            raise ValueError(
                "unexpected number read files to map: %i " % nfiles)

        statement = '''fastq_screen %%(fastq_screen_options)s
                    --outdir %%(outdir)s --conf %%(tempdir)s/fastq_screen.conf
                    %(input_files)s;''' % locals()
        return statement


class Sailfish(Mapper):
    '''run Sailfish to quantify transcript abundance from fastq files'''

    def __init__(self, compress=True, strand="", orient="",
                 threads="", *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.compress = compress

    def mapper(self, infiles, outfile):

        statement = ['''sailfish quant -i %%(index)s''' % locals()]

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        if nfiles == 1:
            input_file = '''-r %s ''' % " ".join(
                ["<(zcat %s)" % x[0] for x in infiles])

        elif nfiles == 2:

            input_file = '''-1 %s -2 %s''' % (
                " ".join(["<(zcat %s)" % x[0] for x in infiles]),
                " ".join(["<(zcat %s)" % x[1] for x in infiles]))

        else:
            # is this the correct error type?
            raise ValueError("incorrect number of input files")

        outdir = os.path.dirname(os.path.abspath(outfile))
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        statement.append('''
        -l %%(sailfish_libtype)s %(input_file)s -o %(outdir)s
        --numBootstraps %%(sailfish_bootstrap)s
        --threads %%(job_threads)s %%(sailfish_options)s;''' % locals())

        statement = " ".join(statement)

        return statement


class Salmon(Mapper):
    '''run Salmon to quantify transcript abundance from fastq files'''

    def __init__(self, compress=True, bias_correct=False,
                 *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.compress = compress

    def mapper(self, infiles, outfile):

        statement = ['''salmon quant -i %%(index)s''' % locals()]

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        if nfiles == 1:
            input_file = '''-r %s ''' % " ".join(
                ["<(zcat %s)" % x[0] for x in infiles])

        elif nfiles == 2:

            input_file = '''-1 %s -2 %s''' % (
                " ".join(["<(zcat %s)" % x[0] for x in infiles]),
                " ".join(["<(zcat %s)" % x[1] for x in infiles]))

        else:
            # is this the correct error type?
            raise ValueError("incorrect number of input files")

        outdir = os.path.dirname(os.path.abspath(outfile))
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        statement.append('''
        -l %%(salmon_libtype)s %(input_file)s -o %(outdir)s
        -k %%(salmon_kmer)s
        --numBootstraps %%(salmon_bootstrap)s
        --threads %%(job_threads)s %%(salmon_options)s;''' % locals())

        statement = " ".join(statement)

        return statement


class Kallisto(Mapper):

    '''run Kallisto to quantify transcript abundance from fastq files
    - set pseudobam to True to output a pseudobam along with the quantification'''

    def __init__(self, pseudobam=False, readable_suffix=False, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

        self.pseudobam = pseudobam
        self.readable_suffix = readable_suffix

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles'''

        tmpdir = os.path.join(self.tmpdir_fastq + "kallisto")

        logfile = outfile + ".log"

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        if nfiles == 1:

            infiles = (" --fragment-length=%(kallisto_fragment_length)s" +
                       " --sd=%(kallisto_fragment_sd)s" +
                       " --single %s" % " ".join([x[0] for x in infiles]))

        elif nfiles == 2:
            infiles = " ".join([" ".join(x) for x in infiles])

        else:
            raise ValueError("incorrect number of input files")

        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        statement = '''
        kallisto quant %%(kallisto_options)s
        -t %%(job_threads)s
        --bootstrap-samples=%%(kallisto_bootstrap)s
        -i %%(index)s -o %(tmpdir)s %(infiles)s''' % locals()

        if self.pseudobam:
            statement += '''
            --pseudobam | samtools view -b -
            > %(outfile)s.bam 2> %(logfile)s;''' % locals()
        else:
            statement += ''' > %(logfile)s &> %(logfile)s ;''' % locals()

        self.tmpdir = tmpdir

        return statement

    def postprocess(self, infiles, outfile):
        '''move outfiles from tmpdir to final location'''

        tmpdir = self.tmpdir

        statement = ('''
        mv -f %(tmpdir)s/abundance.h5 %(outfile)s;
        ''' % locals())

        if self.readable_suffix:

            outfile_readable = outfile + self.readable_suffix
            statement += ('''
            kallisto h5dump -o %(tmpdir)s %(outfile)s;
            mv %(tmpdir)s/abundance.tsv %(outfile_readable)s;
            rm -rf %(tmpdir)s/bs_abundance_*.tsv;''' % locals())

        return statement

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s; rm -rf %s;''' % (
            self.tmpdir_fastq, self.tmpdir)

        return statement


class Counter(Mapper):
    """count number of reads in fastq files."""

    compress = True

    def mapper(self, infiles, outfile):
        '''count number of reads by counting number of lines
        in fastq files.
        '''

        statement = []
        for f in infiles:
            x = " ".join(f)
            statement.append(
                '''zcat %(x)s
                | awk '{n+=1;} END {printf("nreads\\t%%%%i\\n",n/4);}'
                >> %(outfile)s;''' % locals())
        return " ".join(statement)


class SubsetHead(Mapper):
    """subset fastq files by taking the first n sequences"""

    compress = True

    def __init__(self, limit=1000000, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.limit = limit

    def mapper(self, infiles, outfile):
        '''count number of reads by counting number of lines
        in fastq files.
        '''
        limit = self.limit * 4  # 4 lines per fastq entry
        statement = []
        output_prefix = P.snip(outfile, ".subset")
        assert len(infiles) == 1
        infiles = infiles[0]
        if len(infiles) == 1:
            output_filename = output_prefix + ".fastq.gz"
            f = infiles[0]
            statement.append(
                '''zcat %(f)s
                | awk 'NR > %(limit)i {exit} {print}'
                | gzip
                > %(output_filename)s;''' % locals())
        elif len(infiles) > 1:
            for x, f in enumerate(infiles):
                output_filename = output_prefix + ".fastq.%i.gz" % (x + 1)
                statement.append(
                    '''zcat %(f)s
                    | awk 'NR > %(limit)i {exit} {print}'
                    | gzip
                    > %(output_filename)s;''' % locals())
        return " ".join(statement)


class SubsetHeads(Mapper):
    """subset fastq files by taking the first n sequences.

    This differs from SubsetHeads in that list of limits is passed
    along with a list of outfile sentinel file that will be used as a
    file prefix.

    A single file can then be parsed once and subsetted to multiple
    outfiles"""

    compress = True

    def __init__(self, limits=[1000000], *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.limits = limits

    def mapper(self, infiles, outfile):
        '''count number of reads by counting number of lines
        in fastq files.
        '''

        # 4 lines per fastq entry
        limits = [x * 4 for x in sorted(self.limits)]
        output_prefix = P.snip(outfile, ".sentinel")
        assert len(infiles) == 1
        infiles = infiles[0]
        statement = []
        if len(infiles) == 1:
            f = infiles[0]
            awk_cmd = ""
            for n, ix in enumerate(range(0, len(limits))):
                limit = limits[ix]
                output_filename = output_prefix + "_%i.fastq.gz" % n
                awk_cmd += '''{if (NR<%(limit)s) print |
                "gzip > %(output_filename)s"};''' % locals()
            awk_cmd += '{if (NR>%s) {exit}};' % limits[-1]

            statement.append(
                """zcat %(f)s| awk '%(awk_cmd)s';""" % locals())

        elif len(infiles) > 1:
            for x, f in enumerate(infiles):
                awk_cmd = ""
                for n, ix in enumerate(range(0, len(limits))):
                    limit = limits[ix]
                    output_filename = output_prefix + \
                        "_%i.fastq.%i.gz" % (n, x + 1)
                    awk_cmd += '''{if (NR<%(limit)s) print |
                    "gzip > %(output_filename)s"};''' % locals()
                awk_cmd += '{if (NR>%s) {exit}};' % limits[-1]

                statement.append(
                    """zcat %(f)s| awk '%(awk_cmd)s';""" % locals())

        return " ".join(statement)


class SubsetRandom(Mapper):
    """subset fastq files by taking a random n sequences"""

    compress = True

    def __init__(self, limit=1000000, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)
        self.limit = limit

    def mapper(self, infiles, outfile):
        '''count number of reads by counting number of lines
        in fastq files.
        '''
        limit = self.limit * 4  # 4 lines per fastq entry
        statement = []
        output_prefix = P.snip(outfile, ".subset")
        assert len(infiles) == 1
        infiles = infiles[0]

        # check if single or paired end
        if len(infiles) == 1:
            fastq = infiles[0]
            statement = """
            zcat %(fastq1)s |
            paste - - - - |
            sort -R |
            awk -F'\\t'  'NR > %(limit)i {exit} {OFS="\\n";
            print $1,$3,$5,$7 | "gzip > %(output_prefix)s.fastq.gz}"}';
            """ % locals()

        if len(infiles) == 2:
            fastq1, fastq2 = infiles
            statement = """
            paste <(zcat %(fastq1)s) <(zcat %(fastq2)s) |
            paste - - - - |
            sort -R |
            awk -F'\\t'  'NR > %(limit)i {exit} {OFS="\\n";
            print $1,$3,$5,$7 | "gzip > %(output_prefix)s.fastq.1.gz";
            print $2,$4,$6,$8 | "gzip > %(output_prefix)s.fastq.2.gz"}';
            """ % locals()

        return statement


class BWA(Mapper):
    '''Mapper for BWA'''

    def __init__(self,
                 set_nh=False, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

        self.set_nh = set_nh

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with BWA.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using BWA.


        .. note::
            Colour space mapping is not implemented
        '''
        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq + "bwa")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(bwa_aln_options)s"]
        if self.datatype == "solid":
            data_options.append("-c")
            index_prefix = "%(bwa_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            index_prefix = "%(bwa_index_dir)s/%(genome)s"
        else:
            index_prefix = "%(bwa_index_dir)s/%(genome)s"

        data_options = " ".join(data_options)

        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip(os.path.basename(outfile), ".bam")

        if nfiles == 1:
            infiles = " ".join([self.quoteFile(x[0]) for x in infiles])

            statement.append('''
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i
            %(index_prefix)s %(infiles)s
            > %(tmpdir)s/%(track)s.sai 2>>%(outfile)s.bwa.log; checkpoint;
            bwa samse %%(bwa_samse_options)s %%(bwa_index_dir)s/%%(genome)s
            %(tmpdir)s/%(track)s.sai %(infiles)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>>%(outfile)s.bwa.log;
            ''' % locals())

        elif nfiles == 2:
            track1 = track + ".1"
            track2 = track + ".2"
            infiles1 = ",".join([self.quoteFile(x[0]) for x in infiles])
            infiles2 = ",".join([self.quoteFile(x[1]) for x in infiles])

            statement.append('''
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i
            %(index_prefix)s %(infiles1)s
            > %(tmpdir)s/%(track1)s.sai 2>>%(outfile)s.bwa.log; checkpoint;
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i
            %(index_prefix)s %(infiles2)s
            > %(tmpdir)s/%(track2)s.sai 2>>%(outfile)s.bwa.log; checkpoint;
            bwa sampe %%(bwa_sampe_options)s %(index_prefix)s
                      %(tmpdir)s/%(track1)s.sai %(tmpdir)s/%(track2)s.sai
                      %(infiles1)s %(infiles2)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>>%(outfile)s.bwa.log;
            ''' % locals())
        else:
            raise ValueError(
                "unexpected number read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after mapping with BWA.
        This statement removes non-unique reads, sets the NH tag,
        strips sequences and sorts and indexes bam files.

        Parameters
        ----------
        infile: str
            filename of reads file
            can be :term:`fastq`, :term:`sra`, csfasta

        outfile: str
            :term:`bam` filename containing the mapped reads in bam format.

        Returns
        -------
        statement: str
            statement to pass to P.run() to run post processing.
        '''

        # note, this postprocess method is inherited by multiple mappers

        track = P.snip(os.path.basename(outfile), ".bam")
        outf = P.snip(outfile, ".bam")
        tmpdir = self.tmpdir

        strip_cmd, unique_cmd, set_nh_cmd = "", "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique,mapped
            --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence
            --log=%(outfile)s.log''' % locals()

        if self.set_nh:
            set_nh_cmd = '''| cgat bam2bam
            --method=set-nh
            --log=%(outfile)s.log''' % locals()

        statement = '''
                cat %(tmpdir)s/%(track)s.bam
                %(unique_cmd)s
                %(strip_cmd)s
                %(set_nh_cmd)s
                | samtools sort -o %(outfile)s 2>>%(outfile)s.bwa.log;
                samtools index %(outfile)s;''' % locals()

        return statement


class BWAMEM(BWA):

    ''' Mapper for BWA mem algorithm.'''

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with BWA mem.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using BWA.


        .. note::
            Colour space mapping is not implemented
        '''

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq + "bwa")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(bwa_mem_options)s"]
        if self.datatype == "solid":
            data_options.append("-c")
            index_prefix = "%(bwa_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            index_prefix = "%(bwa_index_dir)s/%(genome)s"
        else:
            index_prefix = "%(bwa_index_dir)s/%(genome)s"

        data_options = " ".join(data_options)

        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip(os.path.basename(outfile), ".bam")

        if nfiles == 1:
            infiles = ",".join([self.quoteFile(x[0]) for x in infiles])

            statement.append('''
            bwa mem %%(bwa_mem_options)s -t %%(bwa_threads)i
            %(index_prefix)s %(infiles)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>>%(outfile)s.bwa.log;
            ''' % locals())

        elif nfiles == 2:
            infiles1 = ",".join([self.quoteFile(x[0]) for x in infiles])
            infiles2 = ",".join([self.quoteFile(x[1]) for x in infiles])

            statement.append('''
            bwa mem %%(bwa_mem_options)s -t %%(bwa_threads)i
            %(index_prefix)s %(infiles1)s
            %(infiles2)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>>%(outfile)s.bwa.log;
            ''' % locals())
        else:
            raise ValueError(
                "unexpected number read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)


class Bismark(Mapper):
    '''run bismark to map reads against genome.

    .. note::
       Colour space not implemented

    Arguments
    ---------
    align_stas : bool
        If True, output alignment statistics.
    dedup : bool
        If True, perform a deduplication step.
    read_group_header : bool
        unused?
    '''

    def __init__(self,
                 align_stats=False,
                 dedup=False,
                 read_group_header=False,
                 *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

        self.align_stats = align_stats
        self.dedup = dedup

    def mapper(self, infiles, outfile):

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq + "bwa")
        tmpdir = "."

        bismark_index = "%(bismark_index_dir)s/%(bismark_genome)s"

        tmpdir_fastq = self.tmpdir_fastq
        tmpdir_fastq = "."

        if nfiles == 1:
            infiles = infiles[0][0]
            statement = '''
            bismark %%(bismark_options)s -q --bowtie2
            --output_dir %(tmpdir_fastq)s
            -p %%(bismark_threads)s --bam --phred33-quals %(bismark_index)s
            %(infiles)s;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = infiles[0][0]
            infiles2 = infiles[0][1]

            statement = '''
            bismark %%(bismark_options)s -q --bowtie2
            --output_dir %(tmpdir_fastq)s
            -p %%(bismark_threads)s --bam --non_directional
            --phred33-quals %(bismark_index)s -1 %(infiles1)s -2 %(infiles2)s;
            ''' % locals()

        else:
            raise ValueError(
                "unexpected number read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return statement

    # Postprocessing can identify .sra input from the infiles, then use
    # mapfiles to identify whether the input was single or paired end, then
    # remove reads with low CHH/CHG conversions (i.e more than one unconverted)
    # and rename the bismark outfile to fit ruffus expectations.
    def postprocess(self, infiles, mapfiles, outfile):
        tmpdir_fastq = self.tmpdir_fastq
        tmpdir_fastq = "."
        track = P.snip(os.path.basename(outfile), ".bam")
        infile = infiles[0]
        base = os.path.basename(infile).split(".")[0]
        if infile.endswith(".fastq.gz"):
            statement = '''samtools view -h
            %(tmpdir_fastq)s/%(base)s.fastq.gz_bismark_bt2.bam
            | awk -F" " '$14!~/^XM:Z:[zZhxUu\.]*[HX][zZhxUu\.]*[HX]/ ||
            $1=="@SQ" || $1=="@PG"' | samtools view -b - >
            %%(outdir)s/%(track)s.bam; checkpoint;
            mv %(tmpdir_fastq)s/%(base)s.fastq.gz_bismark_bt2_SE_report.txt
            %%(outdir)s/%(track)s_bismark_bt2_SE_report.txt;''' % locals()
        elif infile.endswith(".fastq.1.gz"):
            statement = '''samtools view -h
            %(tmpdir_fastq)s/%(base)s.fastq.1.gz_bismark_bt2_pe.bam
            | awk -F" " '$14!~/^XM:Z:[zZhxUu\.]*[HX][zZhxUu\.]*[HX]/ ||
            $1=="@SQ" || $1=="@PG"' | samtools view -b - >
            %%(outdir)s/%(track)s.bam; checkpoint;
            mv %(tmpdir_fastq)s/%(base)s.fastq.gz_bismark_bt2_PE_report.txt
            %%(outdir)s/%(track)s_bismark_bt2_SE_report.txt;''' % locals()
        elif infile.endswith(".sra"):
            # this should use Sra module to identify single or paired end
            for mapfile in mapfiles:
                mapfile = os.path.basename(mapfile[0])
                if mapfile.endswith(".fastq.1.gz"):
                    statement = '''samtools view -h
                    %(tmpdir_fastq)s/%(mapfile)s_bismark_bt2_pe.bam|
                    awk -F" " '$14!~/^XM:Z:[zZhxUu\.]*[HX][zZhxUu\.]*[HX]/ ||
                    $1=="@SQ" || $1=="@PG"' | samtools view -b - >
                    %%(outdir)s/%(track)s.bam;checkpoint;
                    mv %(tmpdir_fastq)s/%(mapfile)s_bismark_bt2_PE_report.txt
                    %%(outdir)s/%(track)s_PE_report.txt;''' % locals()
                else:
                    statement = '''samtools view -h
                    %(tmpdir_fastq)s/%(mapfile)s_bismark_bt2.bam|
                    awk -F" " '$14!~/^XM:Z:[zZhxUu\.]*[HX][zZhxUu\.]*[HX]/||
                    $1=="@SQ" || $1=="@PG"' | samtools view -b - >
                    %%(outdir)s/%(track)s.bam;checkpoint;
                    mv %(tmpdir_fastq)s/%(mapfile)s_bismark_bt2_SE_report.txt
                    %%(outdir)s/%(track)s_SE_report.txt;''' % locals()
        else:
            # shouldn't arrive here
            statement = "checkpoint;"

        return statement

    # redefine build method so that mapfiles are
    # also passed to postprocess method
    # this is neccessary for the reasons explained above
    def build(self, infiles, outfile):
        '''run mapper.'''

        cmd_preprocess, mapfiles = self.preprocess(infiles, outfile)
        cmd_mapper = self.mapper(mapfiles, outfile)
        cmd_postprocess = self.postprocess(infiles, mapfiles, outfile)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";")
        assert cmd_mapper.strip().endswith(";")
        if cmd_postprocess:
            assert cmd_postprocess.strip().endswith(";")
        if cmd_clean:
            assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_mapper,
                                          cmd_postprocess,
                                          cmd_clean))

        return statement


class Stampy(BWA):
    '''
    Mapper for Stampy
    '''

    # compress fastq files with gzip
    compress = True

    executable = "stampy.py"

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Stampy.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Stampy.
        '''

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)
        executable = self.executable

        tmpdir = os.path.join(self.tmpdir_fastq + "stampy")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(stampy_align_options)s"]
        if self.datatype == "solid":
            data_options.append("-c")
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s"
        else:
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s"

        track = P.snip(os.path.basename(outfile), ".bam")

        if nfiles == 1:
            infiles = ",".join([self.quoteFile(x[0]) for x in infiles])

            statement.append('''
            %(executable)s -v 3
            -g %%(stampy_index_dir)s/%%(genome)s
            -h %%(stampy_index_dir)s/%%(genome)s
            %%(stampy_options)s
            -M %(infiles)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>%(outfile)s.log;
            ''' % locals())

        elif nfiles == 2:
            track1 = track + ".1"
            track2 = track + ".2"
            infiles1 = ",".join([self.quoteFile(x[0]) for x in infiles])
            infiles2 = ",".join([self.quoteFile(x[1]) for x in infiles])

            statement.append('''
            %(executable)s -v 3 -g %%(stampy_index_dir)s/%%(genome)s
            -h %%(stampy_index_dir)s/%%(genome)s
            %%(stampy_options)s
            -M %(infiles1)s %(infiles2)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam 2>%(outfile)s.log;
            ''' % locals())
        else:
            raise ValueError(
                "unexpected number of read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)


class Butter(BWA):
    '''
    Mapper for Butter
    '''

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Butter.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Butter.
        '''

        if len(infiles) > 1:
            raise ValueError(
                "butter can only operate on one fastq file at a time")

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq + "butter")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip(os.path.basename(outfile), ".bam")

        if nfiles == 1:
            infiles = infiles[0][0]
            track_infile = ".".join(infiles.split(".")[:-1])

            # butter cannot handle compressed fastqs
            # recognises file types by suffix
            track_fastq = os.path.join(tmpdir_fastq, track + ".fastq")

            if infiles.endswith(".gz"):
                statement.append('''
                zcat %(infiles)s > %(track_fastq)s; ''' % locals())
            else:
                statement.append('''
                cat %(infiles)s > %(track_fastq)s; ''' % locals())

            statement.append('''
            butter %%(butter_options)s
            %(track_fastq)s
            %%(butter_index_dir)s/%%(genome)s.fa
            --aln_cores=%%(job_threads)s
            --bam2wig=none
            > %(outfile)s_butter.log;
            ''' % locals())

        elif nfiles == 2:
            raise ValueError(
                "Butter does not support paired end reads ")
        else:
            raise ValueError(
                "unexpected number of read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)

    def cleanup(self, outfile):
        statement = '''rm -rf %s %s;''' % (self.tmpdir_fastq, self.tmpdir)

        return statement


class Tophat(Mapper):
    """
    Mapper for Tophat
    """
    # tophat can map colour space files directly
    preserve_colourspace = True

    # newer versions of tophat can work of compressed files
    compress = True

    executable = "tophat"

    def __init__(self, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Tophat.

        Parameters
        ----------
        infiles: list
            contains 3 filenames
        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            :term:`fasta` filename
            suffix .fa, containing the reference transcriptome

        infiles[2]: str
            junction filename
            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using TopHat.
        '''

        executable = self.executable

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir_tophat = os.path.join(self.tmpdir_fastq + "tophat")
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append("--quals --integer-quals --color")
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join(data_options)

        if nfiles == 1:
            infiles = ",".join([x[0] for x in infiles])
            statement = '''
            %(executable)s --output-dir %(tmpdir_tophat)s
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])

            statement = '''
            %(executable)s --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                    --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()
        elif nfiles == 4:
            # this section works both for paired-ended fastq files
            # in color space mapping (separate quality file)
            # reads1 reads2 qual1 qual2
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])
            infiles3 = ",".join([x[2] for x in infiles])
            infiles4 = ",".join([x[3] for x in infiles])

            statement = '''%(executable)s
                   --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s
                   %(infiles3)s %(infiles4)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        self.tmpdir_tophat = tmpdir_tophat

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after mapping with tophat.
        This statement gzips junctions file, moves log files to output folder,
        indexes using samtools.

        Parameters
        ----------
        infiles: list
            contains three filenames -

        infiles[0]: str
            :term:`fastq` filename
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            :term:`fasta` filename
            suffix .fa, containing the reference transcriptome

        infiles[2]: str
            junction filename
            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
            :term:`bam` filename
            name of the output file from mapping

        Returns
        -------
        statement: str
            statement to pass to P.run() to run post processing.
        '''

        track = P.snip(outfile, ".bam")
        tmpdir_tophat = self.tmpdir_tophat

        statement = '''
        gzip < %(tmpdir_tophat)s/junctions.bed
        > %(track)s.junctions.bed.gz;
        mv %(tmpdir_tophat)s/logs %(outfile)s.logs;
        mv %(tmpdir_tophat)s/accepted_hits.bam %(outfile)s;
        samtools index %(outfile)s;
        ''' % locals()

        return statement


class Tophat2(Tophat):
    '''
    Mapper for tophat2
    '''
    executable = "tophat2"

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Tophat2.

        Parameters
        ----------
        infiles: list
            contains 3 filenames

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            :term:`fasta` filename
            suffix .fa, containing the reference transcriptome

        infiles[2]: str
            junction filename
            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run() to run post processing.
        '''

        # get tophat statement
        statement = Tophat.mapper(self, infiles, outfile)

        # replace tophat options with tophat2 options
        statement = re.sub("%\(tophat_", "%(tophat2_", statement)

        return statement


class Tophat2_fusion(Tophat2):

    executable = "tophat2"

    def mapper(self, infiles, outfile):
        '''build mapping statement for infiles.'''

        # get tophat statement
        statement = Tophat.mapper(self, infiles, outfile)

        # replace tophat options with tophat2 options and add fusion search cmd
        statement = re.sub(re.escape("%(tophat_options)s"),
                           ("%(tophat2_fusion_options)s --fusion-search "
                            "--bowtie1 --no-coverage-search"),
                           statement)

        # replace all other tophat parameters with tophat2_fusion parameters
        statement = re.sub(re.escape("%(tophat_"),
                           "%(tophat2_fusion_", statement)

        return statement

    def postprocess(self, infiles, outfile):
        '''collect output data and postprocess.'''

        track = P.snip(outfile, ".bam")
        tmpdir_tophat = self.tmpdir_tophat

        statement = '''
        gzip < %(tmpdir_tophat)s/junctions.bed
        > %(track)s.junctions.bed.gz;checkpoint;
        mv %(tmpdir_tophat)s/logs %(outfile)s.logs;checkpoint;
        mv %(tmpdir_tophat)s/accepted_hits.bam %(outfile)s;checkpoint;
        mv %(tmpdir_tophat)s/fusions.out %%(fusions)s;checkpoint;
        samtools index %(outfile)s;
        ''' % locals()

        return statement


class TopHat_fusion(Mapper):

    # tophat can map colour space files directly
    preserve_colourspace = False

    def mapper(self, infiles, outfile):
        '''build mapping statement on infiles.
        '''

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir_tophat = os.path.join(self.tmpdir_fastq + "tophat")
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append("--bowtie1 --quals --integer-quals --color")
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join(data_options)

        if nfiles == 1:
            infiles = ",".join([x[0] for x in infiles])
            statement = '''

            tophat2 --output-dir %(tmpdir_tophat)s
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])

            statement = '''

            tophat2 --output-dir %(tmpdir_tophat)s
                    --mate-inner-dist %%(tophat_mate_inner_dist)i
                    --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                  %%(tophat_options)s
                  %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()
        elif nfiles == 4:
            # this section works both for paired-ended fastq files
            # in color space mapping (separate quality file)
            # reads1 reads2 qual1 qual2
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])
            infiles3 = ",".join([x[2] for x in infiles])
            infiles4 = ",".join([x[3] for x in infiles])

            statement = '''

            tophat2 --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s
                   %(infiles3)s %(infiles4)s
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        self.tmpdir_tophat = tmpdir_tophat

        return statement

    def postprocess(self, infiles, outfile):
        '''collect output data and postprocess.'''

        track = P.snip(outfile, "/accepted_hits.sam")
        tmpdir_tophat = self.tmpdir_tophat

        if not os.path.exists('%s' % track):
            os.mkdir('%s' % track)

        # statement = '''
        #    mv -f %(tmpdir_tophat)s/* %(track)s/;
        #    samtools index %(outfile)s;
        #    ''' % locals()
        statement = '''
            mv -f %(tmpdir_tophat)s/* %(track)s/;
            ''' % locals()
        return statement


class Hisat(Mapper):
    '''Mapper for Hisat'''

    # hisat can work of compressed files
    compress = True

    executable = "hisat2"

    def __init__(self, remove_non_unique=False, strip_sequence=False,
                 stranded=False, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

        self.remove_non_unique = remove_non_unique
        self.strip_sequence = strip_sequence
        self.stranded = stranded

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Hisat.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            junction filename
            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Hisat.
        '''

        executable = self.executable

        num_files = [len(x) for x in infiles]
        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")
        nfiles = max(num_files)

        tmpdir_hisat = os.path.join(self.tmpdir_fastq + "hisat")
        tmpdir_fastq = self.tmpdir_fastq
        track = os.path.basename(outfile)

        # add options specific to data type
        index_prefix = "%(hisat_index_dir)s/%(genome)s"

        if self.stranded:
            strandedness = "--rna-strandness %s" % self.stranded
        else:
            strandedness = ""

        if nfiles == 1:
            infiles = ",".join([x[0] for x in infiles])
            statement = '''
            mkdir %(tmpdir_hisat)s;
            %(executable)s
            --threads %%(hisat_threads)i
            %(strandedness)s
            %%(hisat_options)s
            -x %(index_prefix)s
            -U %(infiles)s
            --known-splicesite-infile %%(junctions)s
            > %(tmpdir_hisat)s/%(track)s
            --novel-splicesite-outfile %(outfile)s_novel_junctions
            2>> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])

            statement = '''
            mkdir %(tmpdir_hisat)s;
            %(executable)s
            --threads %%(hisat_threads)i
            %(strandedness)s
            %%(hisat_options)s
            -x %(index_prefix)s
            -1 %(infiles1)s
            -2 %(infiles2)s
            --known-splicesite-infile %%(junctions)s
            > %(tmpdir_hisat)s/%(track)s
            --novel-splicesite-outfile %(outfile)s_novel_junctions
            2>> %(outfile)s.hisat.log;
            ''' % locals()

        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        self.tmpdir_hisat = tmpdir_hisat

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after mapping with hisat.
        This statement strips sequence if specified,
        sorts and indexes with samtools.

        Parameters
        ----------
        infiles: list
            contains two filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            junction filename
            suffix .junctions containing a list of known
            splice junctions.

        outfile: str
        :term:`bam` filename
            name of the output file from mapping

        Returns
        -------
        statement: str
            statement to pass to P.run() to run post processing.
        '''

        track = os.path.basename(outfile)
        tmpdir_hisat = self.tmpdir_hisat

        strip_cmd = ""

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence
            --log=%(outfile)s.log''' % locals()

        statement = '''
        samtools view -uS %(tmpdir_hisat)s/%(track)s
        %(strip_cmd)s
        | samtools sort - -o %(outfile)s 2>>%(outfile)s.hisat.log;
        checkpoint;
        samtools index %(outfile)s;
        checkpoint;
        rm -rf %(tmpdir_hisat)s;
        ''' % locals()

        return statement


class GSNAP(Mapper):
    ''' Mapper for GSNAP'''
    # tophat can map colour space files directly
    preserve_colourspace = True

    # newer versions of tophat can work of compressed files
    compress = True

    executable = "gsnap"

    def __init__(self, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Gsnap.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            filename of type iit containing all known splice sites

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Gsnap.
        '''

        track = P.snip(os.path.basename(outfile), ".bam")

        executable = self.executable

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = self.tmpdir_fastq

        # add options specific to data type
        # index_dir set by environment variable
        index_prefix = "%(gsnap_mapping_genome)s"

        if nfiles == 1:
            if len(infiles) > 1:
                infiles_list = ",".join([x[0] for x in infiles])
                files = "<( zcat %(infiles_list)s)" % locals()
            else:
                individual_infile = infiles[0][0]
                files = "<(zcat %(individual_infile)s)" % locals()

#            statement = '''
#            zcat %(infiles)s
#            | %(executable)s
#                   --nthreads %%(gsnap_worker_threads)i
#                   --format=sam
#                   --db=%(index_prefix)s
#                   %%(gsnap_options)
#                   > %(tmpdir)s/%(track)s.sam
#                   2> %(outfile)s.log;
#            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join([x[0] for x in infiles])
            infiles2 = ",".join([x[1] for x in infiles])

            # patch for compressed files
            if infiles[0][0].endswith(".gz"):
                files = "<( zcat %(infiles1)s ) <( zcat %(infiles2)s )" % \
                        locals()
            else:
                files = "%(infiles1)s %(infiles2)s" % locals()

        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        statement = '''
        %(executable)s
        --nthreads %%(gsnap_worker_threads)i
        --format=sam
        --db=%(index_prefix)s
        %%(gsnap_options)s
        %(files)s
        | samtools view -bS -
        2> %(outfile)s.log ;
        ''' % locals()

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after mapping with Gsnap.

        This statement removes non-unique reads, strips sequences,
        sorts and indexes bam files.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            filename of type iit containing all known splice sites

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''

        track = P.snip(os.path.basename(outfile), ".bam")
        tmpdir = self.tmpdir_fastq

        strip_cmd, unique_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence --log=%(outfile)s.log''' % locals()

        statement = '''
        cat %(tmpdir)s/%(track)s.bam
        %(unique_cmd)s
        %(strip_cmd)s
        | samtools sort -o %(outfile)s 2>>%(outfile)s.log;
        samtools index %(outfile)s;''' % locals()

        return statement


class STAR(Mapper):
    ''' Mapper for STAR'''
    # tophat can map colour space files directly
    preserve_colourspace = True

    # newer versions of tophat can work of compressed files
    compress = True

    executable = "STAR"

    def __init__(self, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Star.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Star.
        '''
        track = P.snip(os.path.basename(outfile), ".bam")

        executable = self.executable

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = self.tmpdir_fastq

        # add options specific to data type
        # index_dir set by environment variable
        index_prefix = "%(genome)s"

        if nfiles == 1:
            infiles = "<( zcat %s )" % " ".join([x[0] for x in infiles])
            statement = '''
            %(executable)s
            --runMode alignReads
            --runThreadN %%(star_threads)i
            --genomeDir %%(star_index_dir)s/%%(star_mapping_genome)s.dir
            --outFileNamePrefix %(tmpdir)s/
            --outStd SAM
            --outSAMunmapped Within
            %%(star_options)s
            --readFilesIn %(infiles)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam
            2> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = " ".join([x[0] for x in infiles])
            infiles2 = " ".join([x[1] for x in infiles])

            if infiles[0][0].endswith(".gz"):
                compress_option = "--readFilesCommand zcat"
            else:
                compress_option = ""

            files = "%(infiles1)s %(infiles2)s" % locals()

            statement = '''
            %(executable)s
            --runMode alignReads
            --runThreadN %%(star_threads)i
            --genomeDir %%(star_index_dir)s/%%(star_mapping_genome)s.dir
            --outFileNamePrefix %(tmpdir)s/
            --outStd SAM
            --outSAMunmapped Within
            %%(star_options)s
            %(compress_option)s
            --readFilesIn %(files)s
            | samtools view -bS -
            > %(tmpdir)s/%(track)s.bam
            2> %(outfile)s.log;
            ''' % locals()

        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after running Star.
        This statement removes non-unique reads, strips sequences,
        sorts and indexes bam files, moves log files to appropriate directory.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''
        track = P.snip(os.path.basename(outfile), ".bam")
        tmpdir = self.tmpdir_fastq

        strip_cmd, unique_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence --log=%(outfile)s.log''' % locals()

        statement = '''
        cp %(tmpdir)s/Log.std.out %(outfile)s.std.log;
        checkpoint;
        cp %(tmpdir)s/Log.final.out %(outfile)s.final.log;
        checkpoint;
        cp %(tmpdir)s/SJ.out.tab %(outfile)s.junctions;
        checkpoint;
        cat %(tmpdir)s/Log.out >> %(outfile)s.log;
        checkpoint;
        cp %(tmpdir)s/Log.progress.out %(outfile)s.progress;
        checkpoint;
        cat %(tmpdir)s/%(track)s.bam
        %(unique_cmd)s
        %(strip_cmd)s
        | samtools sort -o %(outfile)s 2>>%(outfile)s.log;
        samtools index %(outfile)s;''' % locals()

        return statement


class Bowtie(Mapper):
    '''Mapper for Bowtie'''

    # bowtie can map colour space files directly
    preserve_colourspace = True

    executable = "bowtie"

    # prefix to denote index file
    index_option = ""

    # use SAM output in bowtie
    output_option = "--sam"

    def __init__(self, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Bowtie or Bowtie2.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Bowtie.


        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''

        executable = self.executable

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        # transpose files
        infiles = list(zip(*infiles))

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append("--quals --integer-quals --color")
            if nfiles == 2:
                # single end,
                # second file will colors (unpaired data)
                data_options.append(
                    "--quals %s" % ",".join([self.quoteFile(x)
                                             for x in infiles[1]]))
                nfiles -= 1
            elif nfiles == 4:
                raise NotImplementedError()
                data_options.append(
                    "-Q1 %s -Q2 %s" % (",".join(infiles[2], infiles[3])))
                nfiles -= 2
            else:
                raise ValueError("unexpected number of files")
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            data_options.append("-f")
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        index_option = self.index_option
        output_option = self.output_option
        tool_options = self.tool_options
        data_options = " ".join(data_options)
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join([self.quoteFile(x) for x in infiles[0]])
            statement = '''
            %(executable)s --quiet
            --threads %%(bowtie_threads)i
            %(data_options)s
            %(tool_options)s
            %(index_option)s %(index_prefix)s
            %(infiles)s
            %(output_option)s
            2>%(outfile)s_bowtie.log
            | awk -v OFS="\\t" '{sub(/\/[12]$/,"",$1);print}'
            | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam
            1>&2 2>> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = ",".join([self.quoteFile(x) for x in infiles[0]])
            infiles2 = ",".join([self.quoteFile(x) for x in infiles[1]])

            statement = '''
            %(executable)s --quiet
            --threads %%(bowtie_threads)i
            %(data_options)s
            %(tool_options)s
            %(index_option)s %(index_prefix)s
            -1 %(infiles1)s -2 %(infiles2)s
            %(output_option)s
            2>%(outfile)s_bowtie.log
            | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam
            1>&2 2>> %(outfile)s.log;
            ''' % locals()
        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after running Bowtie.
        This statement removes non-unique reads, strips sequences,
        sorts and indexes bam files, moves output files to correct
        directories.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''

        tmpdir_fastq = self.tmpdir_fastq

        unique_cmd, strip_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence --log=%(outfile)s.log''' % locals()

        statement = '''cat %(tmpdir_fastq)s/out.bam
        | cgat bam2bam
        --method=set-nh
        --log=%(outfile)s.log
        %(unique_cmd)s
        %(strip_cmd)s
        | samtools sort -o %(outfile)s;checkpoint;
        samtools index %(outfile)s;
        ''' % locals()

        return statement


class Bowtie2(Bowtie):
    '''Modifies bowtie statements generated with Bowtie mapper for
    Bowtie2'''

    executable = "bowtie2"

    # prefix to denote index file
    index_option = "-x"

    # output option - default is SAM for botwie2
    output_option = ""


class BowtieTranscripts(Mapper):

    '''Mapper for Bowtie against transcripts'''

    # bowtie can map colour space files directly
    preserve_colourspace = True

    compress = True

    executable = "bowtie"

    def __init__(self, *args, **kwargs):
        Mapper.__init__(self, *args, **kwargs)

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Bowtie or Bowtie2
        against the transcriptome.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Bowtie.

        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''
        track = P.snip(outfile, ".bam")

        executable = self.executable

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        # transpose files
        infiles = list(zip(*infiles))

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append("-f -C")
            if nfiles == 2:
                # single end,
                # second file will colors (unpaired data)
                data_options.append("--quals %s" % ",".join(infiles[1]))
                nfiles -= 1
            elif nfiles == 4:
                data_options.append(
                    "-Q1 <( zcat %s ) -Q2 <( zcat %s)" % (
                        ",".join(infiles[2], infiles[3])))
                nfiles -= 2
            else:
                raise ValueError("unexpected number of files")
            index_prefix = "%(prefix)s_cs"
        else:
            index_prefix = "%(prefix)s"

        data_options = " ".join(data_options)
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join(["<(zcat %s)" % x for x in infiles[0]])
            statement = '''
            %(executable)s --quiet --sam
            --threads %%(bowtie_threads)i
            %(data_options)s
            %%(bowtie_options)s
            %(index_prefix)s
            %(infiles)s
            2>%(outfile)s.log
            | awk -v OFS="\\t" '{sub(/\/[12]$/,"",$1);print}'
            | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam
            1>&2 2>> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = ",".join(["<(zcat %s)" % x for x in infiles[0]])
            infiles2 = ",".join(["<(zcat %s)" % x for x in infiles[1]])

            statement = '''
            %(executable)s --quiet --sam
            --threads %%(bowtie_threads)i
            %(data_options)s
            %%(bowtie_options)s
            %(index_prefix)s
            -1 %(infiles1)s -2 %(infiles2)s
            2>%(outfile)s.log
            | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam
            1>&2 2>> %(outfile)s.log;
            ''' % locals()
        else:
            raise ValueError("unexpected number reads to map: %i " % nfiles)

        return statement

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after running Bowtie.
        This statement removes non-unique reads,
        strips sequences, sorts and indexes
        bam files, moves output files to correct directories.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''

        tmpdir_fastq = self.tmpdir_fastq

        strip_cmd, unique_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence --log=%(outfile)s.log''' % locals()

        statement = '''cat %(tmpdir_fastq)s/out.bam
             %(unique_cmd)s
             %(strip_cmd)s
             | samtools sort -o %(outfile)s 2>>%(track)s.bwa.log;checkpoint;
             samtools index %(outfile)s;
             ''' % locals()

        return statement


class BowtieJunctions(BowtieTranscripts):
    '''
    Mapper for bowtie against splice junctions.
    '''

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after running Bowtie and to
        map reads from junction coordinates to genomic coordinates.

        This statement removes non-unique reads,
        strips sequences, sorts and indexes
        bam files, moves output files to correct directories.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''

        tmpdir_fastq = self.tmpdir_fastq

        strip_cmd, unique_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique
            --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence
            --log=%(outfile)s.log''' % locals()

        statement = '''
        cat %(tmpdir_fastq)s/out.bam
        %(unique_cmd)s
        %(strip_cmd)s
        | cgat bam2bam
        --method=set-nh
        --log=%(outfile)s.log
        | cgat rnaseq_junction_bam2bam
        --contigs-tsv-file=%%(contigsfile)s
        --log=%(outfile)s.log
        | samtools sort -o %(outfile)s;
        checkpoint;
        samtools index %(outfile)s;
        ''' % locals()

        return statement


def splitGeneSet(infile):
    ''' split a gtf file by the first column '''

    last = None
    outfile = None
    outprefix = P.snip(infile, ".gtf.gz")

    for line in IOTools.openFile(infile):

        this = line.split("\t")[0]

        if this == last:
            outfile.write(line)

        else:
            last = this
            if outfile is not None:
                outfile.close()

            outfile = IOTools.openFile("%s.%s.gtf.gz" % (outprefix, this), "w")
            outfile.write(line)


class Shortstack(Bowtie):
    '''
    Mapper for Shortstack
    '''

    def mapper(self, infiles, outfile):
        '''
        Build mapping statement from infiles to map with Shortstack.

        Parameters
        ----------
        infiles: list
            contains 1 filename

        infiles[0]: str
            :term:`fastq` filename
            input fastq file containing unmapped reads

        outfile: str
            :term:`bam` filename
            output file to incorporate into the mapping statement

        Returns
        -------
        statement: str
            statement to pass to P.run to map reads using Shortstack.
        '''

        if len(infiles) > 1:
            raise ValueError(
                "Shortstack can only operate on one fastq file at a time")

        num_files = [len(x) for x in infiles]

        if max(num_files) != min(num_files):
            raise ValueError(
                "mixing single and paired-ended data not possible.")

        nfiles = max(num_files)

        tmpdir = os.path.join(self.tmpdir_fastq, "shortstack")
        statement = ["mkdir -p %s;" % tmpdir]
        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip(os.path.basename(outfile), ".shortstack.bam")
        dir_name = os.path.dirname(outfile)

        if nfiles == 1:
            infiles = infiles[0][0]
            track_infile = ".".join(infiles.split(".")[:-1])

            # Shortstack can now handle compressed fastqs
            # have kept code the same but should change to take out
            # redundant code once program is setup and running
            # recognises file types by suffix
            # before you run shortstack at present you need to load module
            # before running the pipeline
            track_fastq = os.path.join(tmpdir_fastq, track + ".fastq")
            if infiles.endswith(".gz"):
                statement.append('''
                zcat %(infiles)s > %(track_fastq)s; ''' % locals())
            else:
                statement.append('''
                cat %(infiles)s > %(track_fastq)s; ''' % locals())

            statement.append('''
            ShortStack %%(shortstack_options)s
            --readfile %(track_fastq)s
            --genomefile %%(shortstack_index_dir)s/%%(genome)s.fa
            --bowtie_cores=%%(job_threads)s
            --outdir %(dir_name)s/%(track)s;
            ''' % locals())

        elif nfiles == 2:
            raise ValueError(
                "Shortstack does not support paired end reads ")
        else:
            raise ValueError(
                "unexpected number of read files to map: %i " % nfiles)

        self.tmpdir = tmpdir

        return " ".join(statement)

    def postprocess(self, infiles, outfile):
        '''
        Builds statement to tidy up after running shortstack.
        This statement removes non-unique reads, strips sequences,
        sorts and indexes bam files, moves output files to correct
        directories.

        Parameters
        ----------
        infiles: list
            contains 2 filenames

        infiles[0]: str
            input file containing unmapped reads
            can be :term:`fastq`, :term:`sra`, csfasta

        infiles[1]: str
            :term:`fasta` file containing reference genome

        outfile: str
            :term:`bam` filename
            output file to incorporate into the statement

        Returns
        -------
        statement: str
            statement to pass to P.run to run post processing.
        '''

        tmpdir_fastq = self.tmpdir_fastq
        track = P.snip(os.path.basename(outfile), ".shortstack.bam")

        unique_cmd, strip_cmd = "", ""

        if self.remove_non_unique:
            unique_cmd = '''| cgat bam2bam
            --method=filter
            --filter-method=unique --log=%(outfile)s.log''' % locals()

        if self.strip_sequence:
            strip_cmd = '''| cgat bam2bam
            --strip-method=all
            --method=strip-sequence --log=%(outfile)s.log''' % locals()

        statement = '''mv shortstack.dir/%(track)s/%(track)s.bam %(outfile)s;
        cat %(outfile)s
        | cgat bam2bam
        --method=set-nh
        --log=%(outfile)s.log
        %(unique_cmd)s
        %(strip_cmd)s
        | samtools sort -o %(outfile)s;
        checkpoint;
        samtools index %(outfile)s;
        ''' % locals()

        return statement

    def cleanup(self, outfile):
        statement = '''rm -rf %s %s;''' % (self.tmpdir_fastq, self.tmpdir)

        return statement
