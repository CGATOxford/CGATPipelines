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


def trackFilters(filtername, bamfile, tabout):
    '''
    Keeps track of which filters are being used on which bam files and
    generates a statement fragment which records this in a log file.
    Echos the name of the filter (filtername) followed by the number of reads
    in the filtered bam file to the table tabout.

    Example Fragment:
    echo unpaired >> K9-13-2_counts.tsv;
    samtools view -c ctmpPGGJro.bam >> K9-13-2_counts.tsv;

    Parameters
    ----------
    filtername: str
        name of filter applied, can be any string
    bamfile: str
        path to filtered bam file
    tabout: str
        path to table to write the output to
    '''
    return """echo %(filtername)s >> %(tabout)s;
              samtools view -c %(bamfile)s.bam >> %(tabout)s; """ % locals()


def appendSamtoolsFilters(statement, inT, tabout, filters, qual, pe):
    '''
    Appends a fragment to an existing command line statement to
    apply filters using samtools
    see https://samtools.github.io/hts-specs/SAMv1.pdf

    Samtools uses the following arguments:
    -F - remove these flags
    -f - keep only these flags

    If the following strings are in the "filters" list these samtools filters
    are applied
    unpaired: -f1, -f2 remove reads which are unpaired or not properly paired
    unmapped: -F4 removes unmapped reads
    secondary: -F 0x100 removed secondary alignments
    lowqual: -q qual removes reads with quality scores < qual

    Example Fragment:
    samtools view -b -q 40 ctmpPGGJro.bam > ctmpnV2rQY.bam;
    rm -f ctmpPGGJro.bam; rm -f ctmpPGGJro;

    The original input file is deleted on the assumption that this is part of
    a list of filters which are applied to a series of temporary files

    Parameters
    ----------
    statement: str
        statement to append to
    inT: str
        path to temporary input bam file - output of previous filtering step
    tabout: str
        path to table to store the number of reads remaining after filtering
    filters: list
       list of filters to apply - this list is searched for
       unmapped, unpaired, secondary, lowqual
    pe: bool
        1 = paired end, 0 = single end
    outfile: str
        path to output file
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

            # filter to a temporary file, remove the original temporary file
            statement += """samtools view -b %(string)s %(inT)s.bam
            > %(outT)s.bam; rm -f %(inT)s.bam; rm -f %(inT)s; """ % locals()
            statement += trackFilters(filt, outT, tabout)
            i += 1

    statement = statement.replace("\n", "")
    return statement, outT


def appendPicardFilters(statement, inT, tabout, filters, pe, outfile):
    '''
    Appends a fragment to an existing command line statement to
    filter bam files using Picard.
    Currently only the MarkDuplicates Picard filter is
    implemented, which removes duplicate reads.

    Example Fragment:

    MarkDuplicates
    INPUT=ctmp87pq2s.bam
    ASSUME_SORTED=true
    REMOVE_DUPLICATES=true
    OUTPUT=ctmphJw7oO.bam
    METRICS_FILE=/dev/null
    VALIDATION_STRINGENCY=SILENT
    2> K9-13-2_filtered_duplicates.log;
    rm -f ctmp87pq2s.bam;
    rm -f ctmp87pq2s;

    The original input file is deleted on the assumption that this is part of
    a list of filters which are applied to a series of temporary files

    Parameters
    ----------
    statement: str
        cmd line statement to append to
    inT: str
        path to temporary input bam file - output of previous filtering step
    tabout: str
        path to table to store the number of reads remaining after filtering
    filters: list
        list of filters to apply - if "duplicates" is in
        this list then the MarkDuplicates fragment will be appended, otherwise
        nothing is appended
    pe: bool
        1 = paired end, 0 = single end
    outfile: str
        path to output file

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
    Appends a fragment to an existing command line statement to
    filters blacklisted regions from a bam file based on a list of bed files
    using either pairToBed (for paired end) or intersect (for single end)
    from bedtools.

    Example Fragment:
    samtools sort -n ctmpnV2rQY.bam -o ctmpL7A2ni.bam;
    rm -f ctmpnV2rQY.bam;
    rm -f ctmpnV2rQY;
    pairToBed -abam ctmpL7A2ni.bam -b chr14.bed -f 0.00000010 -type neither
    > ctmpXbXeki.bam; rm -f ctmpL7A2ni.bam;

    By default a read is removed if it has any overlap with any blacklisted
    region on either end.
    If the read is paired end both halves of the pair are removed if one
    overlaps with a blacklisted region

    The original input file is deleted on the assumption that this is part of
    a list of filters which are applied to a series of temporary files.

    Parameters
    ----------
    statement: str
        cmd line statement to append to
    inT: str
        path to temporary input bam file - output of previous filtering step
    tabout: str
        path to table to store the number of reads remaining after filtering
    blthresh: int
        threshold number of bases of overlap with a blacklisted region
        above which to filter
    pe: bool
        1 = paired end, 0 = single end

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
    Builds a statement which applies various filters to bam files.

    The file is sorted then filters are applied.

    The appendPicardFilters, appendSamtoolsFilters and appendBlacklistFilter
    functions above will add fragments to this statement to
    carry out the filtering steps.

    The filtered bam file is then indexed.
    Read counts after each filtering step are logged in infile_counts.tsv using
    the trackFilters function.

    Example Statement:
    samtools sort K9-13-2.bam -o ctmp87pq2s.bam;
    echo none >> K9-13-2_counts.tsv;
    samtools view -c ctmp87pq2s.bam >> K9-13-2_counts.tsv;
    ...
    samtools sort ctmpXbXeki.bam -o ctmpWm5mq9.bam;
    rm -f ctmpXbXeki.bam;
    rm -f ctmpXbXeki;
    mv ctmpWm5mq9.bam K9-13-2_filtered.bam;
    rm -f ctmpWm5mq9;

    ... represents the statement fragments generated by the filtering functions

    Parameters
    ----------
    infile: str
        path to input file
    outfiles: list
        list of two strings - names of output bam file and output table
    filters: list
        list of filters to apply: can be lowqual, unpaired,
        secondary, unmapped, duplicates
    bedfiles: list
        list of bed files to use as blacklists
    blthresh: int
        threshold base pair overlap with blacklisted regions above which
        to filter
    pe: bool
        1 = paired end, 0 = single end
    strip: bool
        should bam files be stripped after processing, 1 = yes, 0 = no
    qual: int
        minimum mapping quality to keep
    keep_intermediates: bool
        keep temporary files if True

    '''
    bamout, tabout = outfiles
    o = open(tabout, "w")
    o.close()
    inT = infile
    if not os.path.exists("filtered_bams.dir"):
        os.mkdir("filtered_bams.dir")
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
    '''
    Generates a table to ensure that post filtering bam files do not
    contain any of the reads which should have been filtered out.  This table
    is written to a file with the suffix .filteringlog.

    Uses pysam functions is_secondary, is_unmapped, is_proper_pair and
    read.mapping_quality to categorise reads as:
    - proper_pair vs improper_pair
    - high_quality vs low_quality
    - unmapped vs mapped
    - primary vs secondary

    If pipeline_peakcalling is used and filtering has been sucessful then the
    following is expected if these filters are applied:
    "lowqual" filter - low_quality = 0
    "unmapped" filter - unmapped = 0
    "secondary" filter - secondary = 0
    "unpaired" filter - improper_pair = 0

    For paired end, properly paired reads, fragment length is also calculated
    using the pysam template_length function and a second output table is
    generated showing the frequency of each possible fragment length in the
    bam file.  This file will have the suffix .fraglengths.  For unpaired
    data this file is generated but is blank.

    Parameters
    ----------
    infile: str
        input file
    filters:  list
        list of filters to check for
        - can be lowqual, unpaired, secondary, unmapped, duplicates
    qlim: int
        minimum mapping quality to classify as "high quality"
    outfile: str
        path to output file
    '''

    samfile = pysam.AlignmentFile(infile, 'rb')
    sep = '.'
    logfile = sep.join([outfile, 'filteringlog'])

    counter = collections.Counter()
    fragment_length = collections.Counter()
    d = collections.defaultdict(list)

    if "lowqual" not in filters:
        qlim = 0

    counter['primary'] = 0
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
        else:
            counter['primary'] += 1

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

    if "unpaired" in filters and pe == 1:
        outbam = pysam.AlignmentFile("%s.bam" % outfile, "wb",
                                     template=samfile)
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
        outbam.close()
    else:
        outbam = "%s.bam" % outfile
        shutil.copy(infile, outbam)
    out = IOTools.openFile(infile.replace(".bam", ".fraglengths"), "w")
    out.write("frequency\tfrag_length\n")
    for key in fragment_length:
        out.write("%s\t%d\n" % (key, fragment_length[key]))
    out.close()

    filteringReport(counter, logfile)
    samfile.close()


def filteringReport(counter, logfile):
    '''
    Writes a table of the counts contained in the dictionary
    generated in the checkBams function.

    Parameters
    ----------
    counter: dict
        dictionary where keys are types of read (e.g. unpaired)
        and values are the number of reads of this type.
    logfile: str
        path to file to write the output table
    '''
    logfile = open(logfile, "w")
    for c in counter:
        logfile.write("%s\t%s\n" % (c, counter[c]))
    logfile.close()


def estimateInsertSize(infile, outfile, pe, nalignments, m2opts):
    '''
    Predicts fragment size for a bam file and writes it to a table.

    For single end data the MACS2 predictd function is used,
    for paired-end data Bamtools is used.

    In some BAM files the read length (rlen) attribute
    is not set, which causes MACS2 to predict a 0.

    Thus it is computed here from the CIGAR string if rlen is 0.

    Parameters
    ----------
    infile: str
        path to bam file on which to predict fragment size
    outfile: str
        path to output table
    pe: bool
        1 = paired end, 0 = single end
    nalignments: int
        number of lines from bam file to use to predict
        insert size if data is paired end.
    m2opts: str
        string to append when running macs2 predictd
        containing additional options specific to macs2.  Can be an empty
        string for default options, additional options here:
        https://github.com/taoliu/MACS
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
    Generates pseudo bam files by splitting a bam file into two
    equally sized subfiles.  Each read in the input bam is assigned
    to a pseudo bam at random.
    If reads are paired end both reads in the pair are assigned
    to the same bam file.

    Parameters
    ----------
    infile: str
        path to input bam file
    outfiles: list
        list of paths to the two output bam files
    pe: bool
        1 = paired end, 0 = single end
    randomseed: int
        seed to use to generate random numbers
    filters: list
        list of filters previously applied to the bam file.  If this
        list contains the strings 'unpaired' and 'secondary'
        then reads are assumed to be paired end and a check is performed that
        each pseudo bam file contains exactly twice as many reads as read
        names. Anything else in this list is ignored by this function.
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
    '''
    Parses the peak shift estimate file from the estimateInsertSize
    function and returns the fragment size, which is used in the macs2
    postprocessing steps as the "offset" for bed2table
    Parameters
    ----------
    infile: str
        path to input file
    '''

    with IOTools.openFile(infile, "r") as inf:

        header = inf.next().strip().split("\t")
        values = inf.next().strip().split("\t")

        fragment_size_mean_ix = header.index("fragmentsize_mean")

        fragment_size = int(float(values[fragment_size_mean_ix]))

        return fragment_size


def mergeSortIndex(bamfiles, out):
    '''
    Merge bamfiles into a single sorted, indexed outfile.
    Generates and runs a command line statement.

    Example Statement:
    samtools merge ctmpljriXY.bam K9-13-1_filtered.bam K9-13-2_filtered.bam
    K9-13-3_filtered.bam;
    samtools sort ctmpljriXY.bam -o ctmpYH6llm.bam;
    samtools index ctmpYH6llm.bam;
    mv ctmpYH6llm.bam 13_Heart_pooled_filtered.bam;
    mv ctmpYH6llm.bam.bai 13_Heart_pooled_filtered.bam.bai;

    Parameters
    ----------
    bamfiles: list
        list of paths to bam files to merge
    out: str
        path to output file

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
    Sorts and indexes a bam file.
    Generates and runs a command line statement.

    Example Statement:
    samtools sort K9-10-1_filtered.bam -o ctmpvHoczK.bam;
    samtools index ctmpvHoczK.bam;
    mv ctmpvHoczK.bam K9-10-1_filtered.bam;
    mv ctmpvHoczK.bam.bai K9-10-1_filtered.bam.bai

    The input bam file is replaced by the sorted bam file.

    Parameters
    ----------
    bamfile: str
        path to bam file to sort and index

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
    Makes soft links to an existing bam file and its index - used instead
    of copying files.
    Generates and runs a command line statement.

    Parameters:
    currentname: str
        path to original file
    newname: str
        path to link location
    '''
    cwd = os.getcwd()
    os.system("""
    ln -s %(cwd)s/%(currentname)s %(cwd)s/%(newname)s;
    ln -s %(cwd)s/%(currentname)s.bai %(cwd)s/%(newname)s.bai;
    """ % locals())


def makeLink(currentname, newname):
    '''
    Makes a soft link to an existing file- used instead
    of copying files.
    Generates and runs a command line statement.

    Parameters:
    currentname: str
        path to original file
    newname: str
        path to link location
    '''
    cwd = os.getcwd()
    os.system("""
    ln -s %(cwd)s/%(currentname)s %(cwd)s/%(newname)s;
    """ % locals())


class Peakcaller(object):
    '''
    Base class for peakcallers
    Peakcallers call peaks from a BAM file and then post-process peaks
    to generate a uniform output

    Every new peak caller added requires one or more of
    the following functions:
    callPeaks - builds a command line statement to call peaks
    compressOutput - builds a command line statement to compress peakcalling
    output
    postProcessPeaks - builds a command line statement to postprocess peaks
    preparePeaksForIDR - builds a command line statement to prepare IDR input

    Each peak caller should also have a "summarise" function which generates
    a one line summary of the peakcalling results and outputs this to
    a file - these can later be concatenated into a summary table.

    Attributes
    ----------

    threads: int
        number of threads to use for peakcalling
    paired_end: bool
        1 = paired end, 0 = single end
    tool_options: str
        string to append to the cmd statement to run the peakcaller containing
        tool specific options

    '''
    def __init__(self,
                 threads=1,
                 paired_end=True,
                 tool_options=None):
        self.threads = threads
        self.paired_end = paired_end
        self.tool_options = tool_options

    def callPeaks(self, infile, outfile, controlfile):
        '''
        Build command line statement to call peaks.
        Returns an empty string if no callPeaks function is defined for the
        peakcaller used.
        Parameters
        ----------
        infile: str
            path to input bam file
        outfile: str
            path to standard output file (other outputs will also be created
            with the same stem.
        controlfile: path to control (input) bam file
        '''
        return ""

    def compressOutput(self, infile, outfile, contigsfile, controlfile):
        '''
        Build command line statement to compress peakcalling output files.
        Returns an empty string if no compressOutput function is defined for
        the peakcaller used.
        Parameters
        ----------
        infile: str
           path to bam file
        outfile: str
           path to peakcalling output
        contigsfile: str
           path to tab delimited file with the name of each contig in the
           genome in column 0 and contig lengths in column 1. Used by
           bedGraph2bigwig in generating bigwig formatted output files.
        controlfile: str
           path to the control (input) bam file
        '''
        return ""

    def postProcessPeaks(self, infile, outfile, controlfile,
                         insertsizefile):
        '''
        Build command line statement to postprocess peakcalling output.

        Currently none of these outputs are used downstream so any outputs
        of interest to the user of the specific peak caller can be generated.

        Returns an empty string if no postProcess function is defined for
        the peakcaller used

        Parameters
        ----------
        infile: str
           path to bam file
        outfile: str
           path to peakcalling output
        controlfile: str
           path to the control (input) bam file
        insertsizefile: str
           path to table containing insert size data, with columns
           filename, mode, fragmentsize_mean, fragmentsize_std, tagsize
           generated by the estimateInsertSizes function
           this is parsed and used by bedGraphToBigWig
        '''
        return ""

    def preparePeaksForIDR(self, infile, outfile, idr, idrc,
                           idrsuffix, idrcol):
        '''
        Build command line statement to prepare the IDR input.

        IDR requires a sorted file containing the x best peaks from the
        peakcalling step, but the number of peaks to output and the column
        on which to rank depends on the peakcaller.

        Returns an empty string if no preparePeaksForIDR function is defined
        for the peakcaller used

        Parameters
        ----------
        infile: str
            path to bam file
        outfile: str
            path to peakcalling output
        idr: bool
            1 = IDR is enabled 0 = IDR is disabled
        idrc: int
            number of peaks to keep for IDR for this particular peakcaller
        idrsuffix: str
            output file type from this peakcaller to use for IDR
        idrcol: str
            the name of the column on which to rank the peaks for IDR for this
            peakcaller
        '''
        return ""

    def build(self, infile, outfile, contigsfile=None, controlfile=None,
              insertsizef=None, idr=0, idrc=0, idrsuffix=None, idrcol=None):
        '''
        Runs the above functions and uses these to build a complete command
        line statement to run the peakcaller and process its output.

        Parameters
        ----------
        infile: str
           path to bam file
        outfile: str
           path to peakcalling output
        contigsfile: str
           path to tab delimited file with the name of each contig in the
           genome in column 0 and contig lengths in column 1. Used by
           bedGraph2bigwig in generating bigwig formatted output files.
        controlfile: str
           path to the control (input) bam file
        insertsizefile: str
           path to table containing insert size data, with columns
           filename, mode, fragmentsize_mean, fragmentsize_std, tagsize
           generated by the estimateInsertSizes function
           this is parsed and used by bedGraphToBigWig
        idr: bool
            1 = IDR is enabled 0 = IDR is disabled
        idrc: int
            number of peaks to keep for IDR for this particular peakcaller
        idrsuffix: str
            output file type from this peakcaller to use for IDR
        idrcol: int
            the index (0 based) of the column on which to rank the peaks for
            IDR for this peakcaller
        '''

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

    def summarise(self, infile):
        '''
        Function to run after peaks are called and processed to generate a row
        which can later be appended to a summary table for peakcalling for
        all input files.

        Parameters
        ----------
        infile: str
            path to bam file
        '''
        pass


class Macs2Peakcaller(Peakcaller):
    '''
    Peakcaller subclass to call peaks with macs2 and process the macs2 output.

    Attributes
    ----------
    threads: int
        number of threads to use for peakcalling
    paired_end: bool
        1 = paired end, 0 = single end
    tool_options: str
        string to append to the cmd statement with macs2 specific options
    tagsize: int
        if tag size is known it can be added here, otherwise it is calculated
    '''

    def __init__(self,
                 threads=1,
                 paired_end=True,
                 tool_options=None,
                 tagsize=None):
        super(Macs2Peakcaller, self).__init__(threads, paired_end,
                                              tool_options)
        self.tagsize = tagsize

    def callPeaks(self, infile,  outfile, controlfile=None):
        '''
        Build command line statement fragment to call peaks with macs2.

        Example Statement
        macs2 callpeak --format=BAMPE
        --treatment K9-13-2_filtered_pseudo_2.bam --verbose=10
        --name=macs2.dir/K9-13-2_filtered_pseudo_2.macs2
        --qvalue=0.01 --bdg --SPMR --mfold 10 30 --gsize mm
        --broad --broad-cutoff 0.1
        --control IDR_inputs.dir/K9-IN-1_filtered.bam --tsize 75
        >& K9-13-2_filtered_pseudo_2.macs2;
        mv K9-13-2_filtered_pseudo_2.macs2 K9-13-2_filtered_pseudo_2.macs2_log;

        Output files have the same stem but various suffixes.  Details are
        here: https://github.com/taoliu/MACS
        Briefly:
            .macs2_log
             Raw macs2 log file

            .macs2_treat_pileup.bdg
             Bedgraph file of the fragment pileup in the treatment file

            .macs2_control_lambda.brg
             Bedgraph file of the control lambda

            .macs2_peaks.xls
             Tabular file which contains information about called peaks

            .macs2_peaks.broadPeak or .macs2.peaks.narrowPeak
             bed file of peak locations (plus peak summits for narrowPeak)

            .macs2_peaks.gappedPeak
             bed file of narrow and broad peaks

            .macs2_summits.bed
             bed file of summit locations (narrow peaks only)

        The original location for the logging file (suffixed .macs2)
        is overwritten by the output from passing the macs2 xls to
        bed2table.py to give the final required outfile

        Parameters
        ----------
        infile: str
            path to input bam file
        outfile: str
            path to .macs2 output file
        controlfile: str
           path to control (input) bam file
        '''

        if self.tool_options:
            options = [self.tool_options]
        else:
            options = []
        if controlfile:
            options.append("--control %s" % controlfile)

        # tag size is estimated using BamTools
        if self.tagsize is None:
            self.tagsize = BamTools.estimateTagSize(infile)

        options.append("--tsize %i" % self.tagsize)
        options = " ".join(options)

        # Check the paired_end paramter is correct
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
        '''
        Builds a command line statement to compress macs2 outfiles.
        XLS file is compressed using bgzip and tabix
        bedGraph files are compressed using bedGraphToBigWig

        Example Statement:
        bedGraphToBigWig K9-13-2_filtered_pseudo_2.macs2_treat_pileup.bdg
        assembly.dir/contigs.tsv
        K9-13-2_filtered_pseudo_2.macs2_treat_pileup.bw;
        checkpoint;
        rm -rf K9-13-2_filtered_pseudo_2.macs2_treat_pileup.bdg;
        checkpoint;

        bedGraphToBigWig K9-13-2_filtered_pseudo_2.macs2_control_lambda.bdg
        assembly.dir/contigs.tsv
        K9-13-2_filtered_pseudo_2.macs2_control_lambda.bw;
        checkpoint;
        rm -rf K9-13-2_filtered_pseudo_2.macs2_control_lambda.bdg;
        checkpoint;

        grep -v "^$" < K9-13-2_filtered_pseudo_2.macs2_peaks.xls |
        bgzip > K9-13-2_filtered_pseudo_2.macs2_peaks.xls.gz;
        x=$(zgrep "[#|log]" K9-13-2_filtered_pseudo_2.macs2_peaks.xls.gz |
        wc -l);
        tabix -f -b 2 -e 3 -S $x K9-13-2_filtered_pseudo_2.macs2_peaks.xls.gz;
        checkpoint;
        rm -f K9-13-2_filtered_pseudo_2.macs2_peaks.xls;
        checkpoint;

        Output files from the callPeaks step are compressed as follows:
        .macs2_treat_pileup.bdg > .macs2_treat_pileup.bw
        .macs2_control_lambda.brg > .macs2_control_lambda.bw
        .macs2_peaks.xls > .macs2_peaks.xls.gz and .macs2_peaks.xls.gz.tbi

        Parameters
        ----------
         infile: str
            path to input bam file
        outfile: str
            path to .macs2 output file
        contigsfile: str
           path to tab delimited file with the name of each contig in the
           genome in column 0 and contig lengths in column 1. Used by
           bedGraph2bigwig in generating bigwig formatted output files.
        controlfile: str
           path to control (input) bam file
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
            x=$(zgrep "[#|log]" %(outfile)s_%(suffix)s.gz | wc -l);
            tabix -f -b 2 -e 3 -S $x %(outfile)s_%(suffix)s.gz;
             checkpoint; rm -f %(outfile)s_%(suffix)s''' % locals())

        return "; checkpoint ;".join(statement)

    def postProcessPeaks(self, infile, outfile, controlfile, insertsizefile):
        '''
        Generates a command line statement to postprocess MACS 2 results,
        producing further output files which may be of interest.

        The .xls.gz output table from compressOutput is filtered to remove
        comments and column headings then
        bed2table is used to annotate these with various metrics
        including the centre of the peak and
        the number of reads in the sample bam and control bam at the
        peak position.  This output is stored with the suffix .macs

        The .broadPeaks or .narrowPeaks file is processed to output a table,
        suffix .broadpeaks.macs_peaks.bed or .subpeaks.macs_peaks.bed,
        with these metrics plus a "peak height" column.


        Example Statement
        zcat K9-13-2_filtered_pseudo_2.macs2_peaks.xls.gz |
        awk '$1!="chr"' |
        python /ifs/devel/katherineb/cgat/scripts/bed2table.py --counter=peaks
        --bam-file=K9-13-2_filtered_pseudo_2.bam
        --offset=168
        --control-bam-file=K9-IN-1_filtered.bam
        --control-offset=168
        --output-all-fields
        --output-bed-headers=contig,start,end,interval_id,
        -log10\(pvalue\),fold,-log10\(qvalue\),macs_nprobes,macs_peakname
        --log=K9-13-2_filtered_pseudo_2.macs2.log
        > K9-13-2_filtered_pseudo_2.macs;

        cat macs2.dir/K9-13-2_filtered_pseudo_2.macs2_peaks.broadPeak |
        awk '/Chromosome/ {next; } {printf("%s\t%i\t%i\t%i\t%i\n",
        $1,$2,$3,++a,$4)}' |
        python /ifs/devel/katherineb/cgat/scripts/bed2table.py
        --counter=peaks
        --bam-file=K9-13-2_filtered_pseudo_2.bam
        --offset=168
        --control-bam-file=K9-IN-1_filtered.bam
        --control-offset=168
        --output-all-fields
        --output-bed-headers=contig,start,end,interval_id,Height
        --log=K9-13-2_filtered_pseudo_2.macs2.log
        > K9-13-2_filtered_pseudo_2.broadpeaks.macs_peaks.bed;

        Parameters
        ----------
        infile: str
           path to bam file
        outfile: str
           path to peakcalling output
        controlfile: str
           path to the control (input) bam file
        insertsizefile: str
           path to table containing insert size data, with columns
           filename, mode, fragmentsize_mean, fragmentsize_std, tagsize
           generated by the estimateInsertSizes function
           this is parsed and used by bedGraphToBigWig
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
        '''
        Generates a statement fragment which sorts and filters the macs2
        output to provide an input for IDR.
        This involves taking column "idrcol" (currently recommended - column 8
        from the .narrowPeaks or .broadPeaks file), sorting from highest
        to lowest and taking the top n hits (currently recommended 125000).

        Example Statement
        sort -h -r -k8,8 K9-13-2_filtered_pseudo_2.macs2_peaks.broadPeak |
        head -125000 > K9-13-2_filtered_pseudo_2.macs2_IDRpeaks

        Parameters
        ----------
        outfile: str
            path to output file
        idrc: int
            number of top ranked peaks to keep
        idrsuffix: str
            output file type to use for IDR
        idrcol: str
            the index (0 based) of the column on which to rank the peaks for
            IDR
        '''
        statement = ''
        idrout = "%s_IDRpeaks" % outfile
        narrowpeaks = "%s_peaks.%s" % (outfile, idrsuffix)
        col = idrcol
        statement += '''sort -h -r -k%(col)i,%(col)i %(narrowpeaks)s |
        head -%(idrc)s > %(idrout)s''' % locals()
        return statement

    def summarise(self, infile):
        '''
        Parses the MACS2 logfile to extract:
            fragment_size - fragment size
            fragment_treatment_total - number of reads in the treatment bam
            fragment_treatment_filtered - number of reads in the treatment
            after filtering with macs2
            fragment_control_total - number of reads in the control bam
            fragment_control_filtered - number of reads in the control after
            filtering with macs2
            number of peaks - number of peaks passing QC
        This is written to a table with a single row.

        Parameters
        ---------
        infile : str
            path to peakcalling output file (.macs2 file)
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
    '''
    Generates a table containing the pairs of files to compare for IDR
    analysis.
    Various pairs of files are needed:
         Self-consistency analysis
         Pairs of pseudo-bam files generated from the same replicate are
         compared

         Replicate-consistency analysis
         Every possible pair of replicates for each combination of condition
         and treatment are compared

         Pooled-consistency analysis
         Pairs of pseudo bam files generated from pooled replicates for each
         combination of condition and treatment are compared.

    Generates a table with 6 columns:
        file1: path to first input bam file in the pair
        file2: path to the second input bam file in the pair
        IDR_comparison_type: self_consistency, replicate_consistency,
        pooled_consistency
        output_file: path to output table
        Condition: Condition from design table
        Tissue: Tissue from design table

    Parameters
    ----------
    infiles: list
        list of files containing peaks for each replicate, pseudo replicate
        and pooled pseudo replicate
    outilfe: str
        path to output table
    useoracle: bool
        1 = an "oracle peak list" - the peaks from the pooled bam file for each
        combination of condition and tissue will be written to the table as the
        oracle peak list for downstream use
        0 = "None" will be written to the oracle peak list column of the table
        if this is not required
    df: pandas dataframe
        the design table for the experiment, read using the readDesignTable
        function

    '''
    pseudo_reps = []
    pseudo_pooled = []
    notpseudo_reps = []
    notpseudo_pooled = []

    # Categorise files as "pseudo_pooled", "pseudo_reps", "notpseudo_pooled"
    # and "notpseudo_reps
    # pseudo_pooled - peaks from pseudo bam files generated from pooled
    #                 replicates
    # pseudo_reps - peaks from pseudo bam files generated from individual
    #                 replicates
    # notpseudo_reps - peaks from original bam files for individual replicates
    # notpseudo_pooled - peaks from original bam files pooled across replicates

    for f in infiles:
        if "pseudo" in f and "pooled" in f:
            pseudo_pooled.append(f)
        elif "pseudo" in f:
            pseudo_reps.append(f)
        elif "pooled" in f:
            notpseudo_pooled.append(f)
        else:
            notpseudo_reps.append(f)

    # The "oracle peaks" file is the notpseudo_pooled file for each condition
    # and tissue combination
    oracledict = dict()
    # This loop finds the appropriate oracle peak list for the pseudo_pooled
    # pairs and stores this in a dictionary
    cr_pairs = df['Condition'] + "_" + df['Tissue']
    for cr in cr_pairs:
        for npp in notpseudo_pooled:
            npp1 = npp.split("/")[-1]
            if npp1.startswith(cr):
                oracledict[cr] = npp

    # This loop finds the appropriate oracle peak list for the pseudo_reps
    # and pseudo_pooled pairs and stores this in the dictionary

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

    # Generate the table rows for the pseudo_reps peak lists
    # This is achieved by sorting the file names and taking two items at a
    # time from the list - the paired pseudo replicates will be next to
    # each other in this list
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

    # Generates the table rows for the pseudo_pooled peak lists
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

    # segregate the true replicates (notpseudo_reps) into groups
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

    # generate every possible pair of notpseudo_reps
    reppairs = []
    for k in idrrepdict.keys():
        reppairs += list(itertools.combinations(idrrepdict[k], 2))

    # Generate the table rows for the notpseudo_reps
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

    # Write all the table rows to the output file
    out = IOTools.openFile(outfile, "w")
    out.write(
        "file1\tfile2\tIDR_comparison_type\tOracle_Peak_File\tCondition\tTissue\n")
    for p in pairs:
        out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % p)

    out.close()


def buildIDRStatement(infile1, infile2, outfile,
                      sourcec, unsourcec, soft_idr_thresh, idrPARAMS, options,
                      oraclefile=None, test=False):
    '''
    Constructs a command line statement to run the idr software package -
    https://github.com/nboley/idr
    The idr software fails if less than 20 peaks are identified unless the
    --only-merge-peaks parameters is specified, so if "test" is True, this
    parameter is specified and the peaks are counted.  The output of this
    step can then be used to determine if a full analysis should be
    attempted.

    Parameters
    ----------
    infile1: str
        path to first IDR input bam file
    infile2: str
        path to second IDR input bam file
    outfile: str
        path to main IDR output file
    sourcec: str
        statement to source the environment in which to run IDR software
    unsourcec: str
        statement to revert to the original environment
    soft_idr_thresh:
        threshold below which to discard peaks
    idrPARAMS: dict
        dictionary containing parameters for IDR.
        These are:
        idrsuffix:  suffix for the file type to use for IDR for this peakcaller
        for macs2 this is "narrowPeak" or "broadPeak"
        idrcol: 0 based column index to use to rank peaks for IDR, for macs2
        this is 8
        idrcolname: column name of the column referred to in idrcol, for macs2
        this is p.value
        useoracle: should the "oracle" peak list be used instead of the
        standard peak list - 1 yes 0 no
    options: str
        string containing additional options for the idr software
    oraclefile: str
        path to the oracle peak list, if this is requested
    test: bool
        1 = run IDR test to see if > 20 peaks are identified (as otherwise
        the scripts will fail downstream) 0 = run full IDR analysis
    '''
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
    '''
    Builds a summary table for the output from the IDR software.
    This table has the following columns:
    Output_Filename: Name of the main IDR output file
    Replicate_Type: This can be self_consistency, replicate_consistency or
    pooled_consistency depending on the replicate type.
        self_consistency - consistency between pseudo replicates within one
        true replicate
        replicate_consistency - consistency between true replicates
        pooled_consistency - consistency between pseudo replicates within
        pooled replicates for this tissue and condition combination
    Total_Peaks: Number of peaks called in both samples
    Peaks_Passing_IDR: Number of peaks which passed IDR
    Percentage_Peaks_Passing_IDR: (Peaks_Passing_IDR / Total_Peaks) * 100
    Peaks_Failing_IDR: Total_Peaks - Peaks_Passing_IDR
    IDR_Thresholds: Threshold for peaks to pass IDR QC for this replicate type
    IDR_Transformed_Thresholds: -log10 IDR_Threshold
    Input_File_1: IDR input file 1
    Input_File_2: IDR_input file 2
    Condition: Condition from design table
    Tissue: Tissue from design table
    Oracle_Peak_File: List of peaks from pooled replicates (if requested)
    IDR_Successful: Were enough peaks (20) called for IDR to be run -
    True or False
    Conservative_Peak_List: Is this the conservative peak list - the longest
    of the lists of peaks for each condition and tissue combination for the
    replicate_consistency replicates
    Optimal_Peak_List: Is this the optimal peak list - the longest of the
    lists of peaks for each condition and tissue combination for the
    replicate_consistency and pooled_consistency replicates.


    The input files for this step are:
       the filtered IDR output files in bed format
       the IDR summary tables for each IDR test - these have the columns
       "Total_Peaks", "Peaks_Passing_IDR", "Peaks_Failing_IDR",
       "Percentage_Peaks_Failing_IDR" and "IDR_Successful" which correspond to
       the column names in the output file here described above
       a tsv file containing tab delimited pairs of files on which IDR has been
       performed generated using the makePairsForIDR function above.


    Parameters
    ----------
    infiles: tuple
        This tuple should have the layout:
        (((filtered, table), (filtered, table), (filtered, table)), pairs)
        where each (filtered, table) pair is the filtered bed file and summary
        table for an IDR run and "pairs" is the table specifying on which
        pairs of files IDR has been performed.
    outfile: str
        path to output file
    pooledc: float
        the soft threshold used to filter IDR peaks in the pooled_consistency
        replicates
    selfc: float
        the soft threshold used to filter IDR peaks in the self_consistency
        replicates
    repc: float
        the soft threshold used to filter IDR peaks in the
        replicate_consistency replicates
    '''
    tables = [i[1] for i in infiles[:-1]]
    pairs = infiles[-1]

    # read all the IDR summary tables into a single file
    alltab = pd.DataFrame()
    for table in tables:
        tab = pd.read_csv(table, sep="\t")
        tab['Output_Filename'] = table.split("/")[-1]
        alltab = alltab.append(tab)

    # read the table containing the pairs of files on which IDR
    # has been performed
    pairtab = pd.read_csv(pairs, sep="\t", names=['Input_File_1',
                                                  'Input_File_2',
                                                  'Replicate_Type',
                                                  'Oracle_Peak_File',
                                                  'Condition',
                                                  'Tissue'])

    # Generate a temporary column to merge the two tables corresponding to
    # the "Output_Filename" column in the summary table
    pairstrings = []
    for p in pairtab.index.values:
        p = list(pairtab.ix[p])
        p1 = P.snip(p[0].split("/")[-1])
        p2 = P.snip(p[1].split("/")[-1])
        pairstring = "%s_v_%s_table.tsv" % (p1, p2)
        pairstrings.append(pairstring)

    # Tidy up the column names etc in the output tables
    pairtab['pairstring'] = pairstrings
    alltab = alltab.merge(pairtab, left_on='Output_Filename',
                          right_on='pairstring')
    alltab = alltab.drop('pairstring', 1)

    # Find the appropriate thresholds for each replicate type and store in a
    # list
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

    # find the "Conservative_Peak_List" for each experiment
    alltab['Conservative_Peak_List'] = False
    replicates = alltab[alltab['Replicate_Type'] == 'replicate_consistency']
    for exp in set(replicates['Experiment'].values):
        subtab = replicates[replicates['Experiment'] == exp]
        m = max(subtab['Peaks_Passing_IDR'])
        val = subtab[subtab['Peaks_Passing_IDR'] == m].index.values
        alltab.ix[val, 'Conservative_Peak_List'] = True

    # find the "Optimal_Peak_List" for each experiment
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
    '''
    Generates a table containing various quality control statistics for the IDR
    analysis.
    The input file is the output table from the summariseIDR function
    above.
    The output table has the following columns:
    Experiment: Tissue and condition
    Ratio_Type: Name of this summary statistic
    Ratio_Formula: How the statistic was calculated
    Ratio: Value of summary statistic
    Comparison_File_1: First file used to generate this statistic
    Comparison_File_2: Second file used to generate this statistic.

    The following statistics are calculated:
    self-consistency ratio - this is calculated as N1/N2 for every possible
    pair of true replicates within a condition and tissue combination, where
    N1 and N2 are the number of peaks passing IDR for the self-consistency
    test (comparing two pseudo replicates) for replicate 1 and replicate 2
    respectively.
    
    rescue ratio - this is calculated for each combination of tissue and
    condition as max(Np, Nt) / min(Np, Nt) where
    Np is the number of peaks passing the pooled consistency IDR (comparing
    pseudo replicates generated from the pooled bam file for the experiment)
    and Nt is the length of the longest peak list generated from the
    replicate_consistency tests of this experiment.

    Parameters
    ----------
    infile: str
        path to input table, which should be the output from summariseIDR
        above
    outfile: str
        path to output table
    '''
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
    '''
    Runs the R package ChIPQC to plot and record various quality control
    statistics on peaks.
    
    Parameters
    ----------
    infile: str
       path to a design table formatted for the ChIPQC package as specified
       here  -
       https://www.bioconductor.org/packages/devel/bioc/
       vignettes/ChIPQC/inst/doc/ChIPQC.pdf
    outfile: str
       path to main output file
    rdir: str
       path to directory in which to place the output files
    '''
    
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


# Pipeline Specific Functions


def readDesignTable(infile, poolinputs):
    '''
    This function reads a design table named "design.tsv"  and generates
    objects to be used to match peaks called in samples to the appropriate
    inputs prior to IDR analysis.

    These objects are:

    1. A dictionary, inputD, linking each input file and each of the various
    subfiles required for IDR to the appropriate input,
    as specified in the design table

    2. A pandas dataframe, df, containing the information from the
    design table

    This dictionary includes the filenames that the IDR output files will be
    generated downstream in pipeline_peakcalling, so if these steps are run
    outside of the pipeline they need to be named according to the same
    convention.

    If IDR is not requsted, these dictionary values will still exist but will
    not be used.

    poolinputs has three possible options:

    none
    Use the input files per replicate as specified in the design file
    Input files will still also be pooled if IDR is specified as IDR requires
    BAM files representing pooled replicates as well as BAM files for each
    replicate

    all
    Pool all input files and use this single pooled input BAM file as the input
    for any peakcalling.  Used when input is low depth or when only a single
    input or replicates of a single input are available.

    condition
    pool the input files within each combination of conditions and tissues

    Parameters
    ----------
    infile: path to design file
    poolinputs: can be "none", "all" or "condition" as specified above
    '''

    # read the design table
    df = pd.read_csv(infile, sep="\t")
    pairs = zip(df['bamReads'], df['bamControl'])
    CHIPBAMS = list(set(df['bamReads'].values))

    # if poolinputs is none, the inputs are as specified in the design table
    if poolinputs == "none":
        inputD = dict(pairs)
        conditions = df['Condition'].values
        tissues = df['Tissue'].values
        i = 0
        # preempts what the pooled input files for IDR analysis will be
        # named
        for C in CHIPBAMS:
            cond = conditions[i]
            tissue = tissues[i]
            inputD["%s_%s.bam" % (cond, tissue)] = "%s_%s_pooled.bam" % (
                cond, tissue)
            i += 1

    elif poolinputs == "all":
        # if poolinputs is all, all inputs will be pooled and used any time
        # an input is needed
        inputD = dict()
        for C in CHIPBAMS:
            inputD[C] = "pooled_all.bam"

    elif poolinputs == "condition":
        inputD = dict()
        conditions = df['Condition'].values
        tissues = df['Tissue'].values
        i = 0
        # preempts the name of all input files for all combinations of
        # files while will be generated for the IDR
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


def readTable(tabfile):
    '''
    Used by the pipeline_peakcalling.py pipeline
    to read the "peakcalling_bams_and_inputs.tsv" file back
    into memory.
    Parameters
    ----------
    tabfile: str
        path to peakcalling_bams_and_inputs.tsv table
    '''
    df = pd.read_csv(tabfile, sep="\t")
    chips = df['ChipBam'].values
    inputs = df['InputBam'].values
    pairs = zip(chips, inputs)
    D = dict(pairs)
    return D
