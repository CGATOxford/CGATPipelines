"""
PipelineMappingQC.py - Tasks for QC'ing mapping
===============================================

Reference
---------

"""

import CGAT.Experiment as E
import os
import re
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools
import CGATPipelines.Pipeline as P
import pandas as pd

PICARD_MEMORY = "5G"


def getNumReadsFromReadsFile(infile):
    '''get number of reads from a .nreads file.'''
    with IOTools.openFile(infile) as inf:
        line = inf.readline()
        if not line.startswith("nreads"):
            raise ValueError(
                "parsing error in file '%s': "
                "expected first line to start with 'nreads'")
        nreads = int(line[:-1].split("\t")[1])
    return nreads


def buildPicardInsertSizeStats(infile, outfile, genome_file):
    '''run Picard:CollectInsertSizeMetrics

    Collect insert size statistics.

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.
    '''

    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    statement = '''CollectInsertSizeMetrics
    INPUT=%(infile)s
    REFERENCE_SEQUENCE=%(genome_file)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s'''

    P.run()


def buildPicardAlignmentStats(infile, outfile, genome_file):
    '''run picard:CollectMultipleMetrics

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.
    '''

    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    # Picard seems to have problem if quality information is missing
    # or there is no sequence/quality information within the bam file.
    # Thus, add it explicitly.
    statement = '''cat %(infile)s
    | cgat bam2bam -v 0
    --method=set-sequence --output-sam
    | CollectMultipleMetrics
    INPUT=/dev/stdin
    REFERENCE_SEQUENCE=%(genome_file)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s'''

    P.run()


def buildPicardDuplicationStats(infile, outfile):
    '''run picard:MarkDuplicates

    Record duplicate metrics using Picard, the marked records
    are discarded.

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    '''

    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    # currently, MarkDuplicates cannot handle split alignments from gsnap
    # these can be identified by the custom XT tag.
    if ".gsnap.bam" in infile:
        tmpf = P.getTempFile(".")
        tmpfile_name = tmpf.name
        statement = '''samtools view -h %(infile)s
        | awk "!/\\tXT:/"
        | samtools view /dev/stdin -S -b > %(tmpfile_name)s;
        ''' % locals()
        data_source = tmpfile_name
    else:
        statement = ""
        data_source = infile

    os.environ["CGAT_JAVA_OPTS"] = "-Xmx%s -XX:+UseParNewGC\
                                    -XX:+UseConcMarkSweepGC" % (PICARD_MEMORY)

    statement += '''MarkDuplicates
    INPUT=%(data_source)s
    ASSUME_SORTED=true
    METRICS_FILE=%(outfile)s
    OUTPUT=/dev/null
    VALIDATION_STRINGENCY=SILENT
    '''
    P.run()

    os.unsetenv("CGAT_JAVA_OPTS")

    if ".gsnap.bam" in infile:
        os.unlink(tmpfile_name)


def buildPicardDuplicateStats(infile, outfile):
    '''run picard:MarkDuplicates

    Record duplicate metrics using Picard and keep the dedupped .bam
    file.

    Pair duplication is properly handled, including inter-chromosomal
    cases. SE data is also handled.  These stats also contain a
    histogram that estimates the return from additional sequecing.  No
    marked bam files are retained (/dev/null...)  Note that picards
    counts reads but they are in fact alignments.

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.

    '''
    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    os.environ["CGAT_JAVA_OPTS"] = "-Xmx%s -XX:+UseParNewGC\
                                    -XX:+UseConcMarkSweepGC" % (PICARD_MEMORY)
    statement = '''MarkDuplicates
    INPUT=%(infile)s
    ASSUME_SORTED=true
    METRICS_FILE=%(outfile)s.duplicate_metrics
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT;
    '''
    statement += '''samtools index %(outfile)s ;'''
    P.run()
    os.unsetenv("CGAT_JAVA_OPTS")


def buildPicardCoverageStats(infile, outfile, baits, regions):
    '''run picard:CalculateHSMetrics

    Generate coverage statistics for regions of interest from a bed
    file using Picard.

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    baits : :term:`bed` formatted file of bait regions
    regions : :term:`bed` formatted file of target regions

    '''

    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    statement = '''CalculateHsMetrics
    BAIT_INTERVALS=%(baits)s
    TARGET_INTERVALS=%(regions)s
    INPUT=%(infile)s
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=LENIENT''' % locals()
    P.run()


def buildPicardGCStats(infile, outfile, genome_file):
    """picard:CollectGCBiasMetrics

    Collect GC bias metrics.

    Arguments
    ---------
    infile : string
        Input filename in :term:`BAM` format.
    outfile : string
        Output filename with picard output.
    genome_file : string
        Filename with genomic sequence.

    """

    job_memory = PICARD_MEMORY
    job_threads = 3

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    statement = '''CollectGcBiasMetrics
    INPUT=%(infile)s
    REFERENCE_SEQUENCE=%(genome_file)s
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    CHART_OUTPUT=%(outfile)s.pdf
    SUMMARY_OUTPUT=%(outfile)s.summary
    >& %(outfile)s'''

    P.run()


def loadPicardMetrics(infiles, outfile, suffix,
                      pipeline_suffix=".picard_stats",
                      tablename=None):
    '''load picard metrics.

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``toTable(outfile) + "_" +
        suffix``.

    '''

    if not tablename:
        tablename = "%s_%s" % (P.toTable(outfile), suffix)

    outf = P.getTempFile(".")

    filenames = ["%s.%s" % (x, suffix) for x in infiles]

    first = True
    for filename in filenames:
        track = P.snip(os.path.basename(filename), "%s.%s" %
                       (pipeline_suffix, suffix))

        if not os.path.exists(filename):
            E.warn("File %s missing" % filename)
            continue

        lines = IOTools.openFile(filename, "r").readlines()

        # extract metrics part
        rx_start = re.compile("## METRICS CLASS")
        for n, line in enumerate(lines):
            if rx_start.search(line):
                lines = lines[n + 1:]
                break

        for n, line in enumerate(lines):
            if not line.strip():
                lines = lines[:n]
                break

        if len(lines) == 0:
            E.warn("no lines in %s: %s" % (track, filename))
            continue

        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
            fields = lines[0][:-1].split("\t")
        else:
            f = lines[0][:-1].split("\t")
            if f != fields:
                raise ValueError(
                    "file %s has different fields: expected %s, got %s" %
                    (filename, fields, f))

        first = False
        for i in range(1, len(lines)):
            outf.write("%s\t%s" % (track, lines[i]))

    outf.close()

    P.load(outf.name,
           outfile,
           tablename=tablename,
           options="--add-index=track --allow-empty-file")

    os.unlink(outf.name)


def loadPicardHistogram(infiles, outfile, suffix, column,
                        pipeline_suffix=".picard_stats", tablename=False):
    '''extract a histogram from a picard output file and load
    it into database.

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile.
    suffix : string
        Suffix to append to table name.
    column : string
        Column name to take from the histogram.
    pipeline_suffix : string
        Suffix to remove from track name.
    tablename : string
        Tablename to use. If unset, the table name will be derived
        from `outfile` and suffix as ``toTable(outfile) + "_" +
        suffix``.
    '''

    if not tablename:
        tablename = "%s_%s" % (P.toTable(outfile), suffix)
        tablename = tablename.replace("_metrics", "_histogram")

    # some files might be missing
    xfiles = [x for x in infiles if os.path.exists("%s.%s" % (x, suffix))]

    if len(xfiles) == 0:
        E.warn("no files for %s" % tablename)
        return

    header = ",".join([P.snip(os.path.basename(x), pipeline_suffix)
                       for x in xfiles])
    filenames = " ".join(["%s.%s" % (x, suffix) for x in xfiles])

    # there might be a variable number of columns in the tables
    # only take the first ignoring the rest

    load_statement = P.build_load_statement(
        tablename,
        options="--add-index=track "
        " --header-names=%s,%s"
        " --allow-empty-file"
        " --replace-header" % (column, header))

    statement = """cgat combine_tables
    --regex-start="## HISTOGRAM"
    --missing-value=0
    --take=2
    %(filenames)s
    | %(load_statement)s
    >> %(outfile)s
    """

    P.run()


def loadPicardAlignmentStats(infiles, outfile):
    '''load all output from Picard's CollectMultipleMetrics into database.

    Loads tables into database with prefix derived from outfile:
       * [outfile]_alignment_summary_metric
       * [outfile]_insert_size_metrics
       * [outfile]_quality_by_cycle_metrics
       * [outfile]_quality_distribution_metrics
       * [outfile]_insert_size_metrics

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.

    '''

    loadPicardMetrics(infiles, outfile, "alignment_summary_metrics")

    # insert size metrics only available for paired-ended data
    loadPicardMetrics(infiles, outfile, "insert_size_metrics")

    histograms = (("quality_by_cycle_metrics", "cycle"),
                  ("quality_distribution_metrics", "quality"),
                  ("insert_size_metrics", "insert_size"))

    for suffix, column in histograms:
        loadPicardHistogram(infiles, outfile, suffix, column)


def loadPicardDuplicationStats(infiles, outfiles):
    '''load picard duplicate filtering stats into database.

    Loads two tables into the database
       * picard_duplication_metrics
       * picard_complexity_histogram

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''
    # SNS: added to enable naming consistency

    outfile_metrics, outfile_histogram = outfiles

    suffix = "picard_duplication_metrics"

    # the loading functions expect "infile_name.pipeline_suffix" as the infile
    # names.
    infile_names = [x[:-len("." + suffix)] for x in infiles]

    loadPicardMetrics(infile_names, outfile_metrics, suffix, "",
                      tablename="picard_duplication_metrics")

    infiles_with_histograms = []

    # The complexity histogram is only present for PE data, so we must check
    # because by design the pipeline does not track endedness
    for infile in infile_names:
        with_hist = False
        with open(".".join([infile, suffix]), "r") as open_infile:
            for line in open_infile:
                if line.startswith("## HISTOGRAM"):
                    infiles_with_histograms.append(infile)
                    break

    if len(infiles_with_histograms) > 0:
        loadPicardHistogram(infiles_with_histograms,
                            outfile_histogram,
                            suffix,
                            "coverage_multiple",
                            "",
                            tablename="picard_complexity_histogram")
    else:
        with open(outfile_histogram, "w") as ofh:
            ofh.write("No histograms detected, no data loaded.")


def loadPicardDuplicateStats(infiles, outfile, pipeline_suffix=".bam"):
    '''load picard duplicate filtering stats.

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    pipeline_suffix : string
        Suffix appended to pipeline output file, will be removed to
        define track.
    '''

    loadPicardMetrics(
        infiles, outfile, "duplicate_metrics",
        pipeline_suffix=pipeline_suffix)
    loadPicardHistogram(infiles,
                        outfile,
                        "duplicate_metrics",
                        "duplicates",
                        pipeline_suffix=pipeline_suffix)


def loadPicardCoverageStats(infiles, outfile):
    '''import coverage statistics into database.

    Arguments
    ---------
    infiles : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    outf = P.getTempFile(".")
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".cov")
        lines = [x for x in open(f, "r").readlines()
                 if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
    outf.close()
    P.load(outf.name,
           outfile,
           options="--ignore-empty --add-index=track")
    os.unlink(outf.name)


def buildBAMStats(infile, outfile):
    '''Count number of reads mapped, duplicates, etc. using
    bam2stats.py

    Arguments
    ---------
    infile : string
        Input file in :term:`BAM` format
    outfile : string
        Output file in :term:`tsv` format.

    '''

    statement = '''cgat bam2stats
    --force-output
    --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s'''
    P.run()


def loadBAMStats(infiles, outfile):
    '''load output of :func:`buildBAMStats` into database.

    Arguments
    ---------
    infiles : string
        Input files, output from :func:`buildBAMStats`.
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    header = ",".join([P.snip(os.path.basename(x), ".readstats")
                       for x in infiles])
    filenames = " ".join(["<( cut -f 1,2 < %s)" % x for x in infiles])
    tablename = P.toTable(outfile)

    load_statement = P.build_load_statement(
        tablename,
        options="--add-index=track "
        " --allow-empty-file")

    E.info("loading bam stats - summary")
    statement = """cgat combine_tables
    --header-names=%(header)s
    --missing-value=0
    --ignore-empty
    %(filenames)s
    | perl -p -e "s/bin/track/"
    | cgat table2table --transpose
    | %(load_statement)s
    > %(outfile)s"""
    P.run()

    for suffix in ("nm", "nh"):
        E.info("loading bam stats - %s" % suffix)
        filenames = " ".join(["%s.%s" % (x, suffix) for x in infiles])

        load_statement = P.build_load_statement(
            "%s_%s" % (tablename, suffix),
            options="--allow-empty-file")

        statement = """cgat combine_tables
        --header-names=%(header)s
        --skip-titles
        --missing-value=0
        --ignore-empty
        %(filenames)s
        | perl -p -e "s/bin/%(suffix)s/"
        | %(load_statement)s
        >> %(outfile)s """
        P.run()

    # load mapping qualities, there are two columns per row
    # 'all_reads' and 'filtered_reads'
    # Here, only filtered_reads are used (--take=3)
    for suffix in ("mapq",):
        E.info("loading bam stats - %s" % suffix)
        filenames = " ".join(["%s.%s" % (x, suffix) for x in infiles])

        load_statement = P.build_load_statement(
            "%s_%s" % (tablename, suffix),
            options=" --allow-empty-file")

        statement = """cgat combine_tables
        --header-names=%(header)s
        --skip-titles
        --missing-value=0
        --ignore-empty
        --take=3
        %(filenames)s
        | perl -p -e "s/bin/%(suffix)s/"
        | %(load_statement)s
        >> %(outfile)s """
        P.run()


def buildPicardRnaSeqMetrics(infiles, strand, outfile):
    '''run picard:RNASeqMetrics



    Arguments
    ---------
    infiles : string
        Input filename in :term:`BAM` format.
        Genome file in refflat format
            (http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat)
    outfile : string
        Output filename with picard output.

    '''
    job_memory = PICARD_MEMORY
    job_threads = 3
    infile, genome = infiles

    if BamTools.getNumReads(infile) == 0:
        E.warn("no reads in %s - no metrics" % infile)
        P.touch(outfile)
        return

    os.environ["CGAT_JAVA_OPTS"] = "-Xmx%s -XX:+UseParNewGC\
                                    -XX:+UseConcMarkSweepGC" % (PICARD_MEMORY)
    statement = '''CollectRnaSeqMetrics
    REF_FLAT=%(genome)s
    INPUT=%(infile)s
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    STRAND=%(strand)s
    VALIDATION_STRINGENCY=SILENT
    '''
    P.run()
    os.unsetenv("CGAT_JAVA_OPTS")


def loadPicardRnaSeqMetrics(infiles, outfiles):
    '''load picard rna stats into database.

    Loads tables into the database
       * picard_rna_metrics
       * picard_rna_histogram

    Arguments
    ---------
    infile : string
        Filenames of files with picard metric information. Each file
        corresponds to a different track.
    outfiles : string
        Logfile. The table names will be derived from `outfile`.
    '''

    outfile_metrics, outfile_histogram = outfiles

    suffix = "picard_rna_metrics"

    # the loading functions expect "infile_name.pipeline_suffix" as the infile
    # names.
    infile_names = [x[:-len("." + suffix)] for x in infiles]

    loadPicardMetrics(infile_names, outfile_metrics, suffix, "",
                      tablename="picard_rna_metrics")

    infiles_with_histograms = []

    # Checking if histogram is present (?is this necessary)
    for infile in infile_names:
        with_hist = False
        with open(".".join([infile, suffix]), "r") as open_infile:
            for line in open_infile:
                if line.startswith("## HISTOGRAM"):
                    infiles_with_histograms.append(infile)
                    break

    if len(infiles_with_histograms) > 0:
        loadPicardHistogram(infiles_with_histograms,
                            outfile_histogram,
                            suffix,
                            "coverage_multiple",
                            "",
                            tablename="picard_rna_histogram")
    else:
        with open(outfile_histogram, "w") as ofh:
            ofh.write("No histograms detected, no data loaded.")


def loadIdxstats(infiles, outfile):
    '''take list of file paths to samtools idxstats output files
    and merge to create single dataframe containing mapped reads per
    contig for each track. This dataframe is then loaded into
    database.

    Loads tables into the database
        * idxstats_reads_per_chromosome

    Arguments
    ---------
    infiles : list
        list where each element is a string of the filename containing samtools
        idxstats output. Filename format is expected to be 'sample.idxstats'
    outfile : string
        Logfile. The table name will be derived from `outfile`.
    '''

    outf = P.getTempFile(".")
    dfs = []
    for f in infiles:
        track = P.snip(f, ".idxstats").split('/')[-1]

        if not os.path.exists(f):
            E.warn("File %s missing" % f)
            continue

        # reformat idx stats
        df = pd.read_csv(f, sep='\t', header=None)
        df.columns = ['region', 'length', 'mapped', 'unmapped']

        # calc total reads mapped & unmappedpep
        total_reads = df.unmapped.sum() + df.mapped.sum()
        total_mapped_reads = df.mapped.sum()

        reformatted_df = pd.DataFrame([['total_mapped_reads', total_mapped_reads],
                                      ['total_reads', total_reads],
                                      ['track', track]], columns=(['region', 'mapped']))

        # reformat the df
        df = df.append(reformatted_df, ignore_index=True)
        df.set_index('region', inplace=True)
        df1 = df[['mapped']].T
        dfs.append(df1)

    # merge dataframes into single table
    master_df = pd.concat(dfs)
    master_df.drop('*', axis=1, inplace=True)
    master_df.to_csv(outf, sep='\t', index=False)
    outf.close()

    P.load(outf.name,
           outfile,
           options="--ignore-empty --add-index=track")
    os.unlink(outf.name)
