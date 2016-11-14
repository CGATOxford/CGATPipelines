"""PipelineReadqc.py - Tasks for QC'ing short read data sets
============================================================

The majority of the functions in this module are for running
and processing the information from the fastqc_ tool.

Reference
---------

"""

import os
import re
import glob
import collections
from six import StringIO
import pandas as pd
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.CSV2DB as CSV2DB


def FastqcSectionIterator(infile):
    """iterate over FASTQC output file and yield each section.

    Sections in FASTQC output start with `>>`` and end with
    ``>>END_MODULE``.

    Yields
    ------
    name : string
        Section name
    status : string
        Section status
    header : string
        Section header
    data : list
        Lines within section

    Arguments
    ---------
    infile : iterator
        Iterator over contents of Fastqc output.

    """
    data = []
    name, status, header, data = None, None, None, None
    for line in infile:
        if line.startswith(">>END_MODULE"):
            yield name, status, header, data
        elif line.startswith(">>"):
            name, status = line[2:-1].split("\t")
            data = []
        elif line.startswith("#"):
            header = "\t".join([x for x in line[1:-1].split("\t") if x != ""])
        else:
            data.append(
                "\t".join([x for x in line[:-1].split("\t") if x != ""]))


def collectFastQCSections(infiles, section, datadir):
    '''iterate over all fastqc files and extract a particular section.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    section : string
        Section name to extract
    datadir : string
        Location of actual Fastqc output to be parsed.

    Returns
    -------
    results : list
        List of tuples, one tuple per input file. Each tuple contains
        track, status, header and data of `section`.

    '''
    results = []
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".fastqc")
        filename = os.path.join(datadir, track + "*_fastqc", "fastqc_data.txt")
        for fn in glob.glob(filename):
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                if name == section:
                    results.append((track, status, header, data))
    return results


def loadFastqc(filename,
               backend="sqlite",
               database="csvdb",
               host="",
               username="",
               password="",
               port=3306):
    '''load FASTQC statistics into database.

    Each section will be uploaded to its own table.

    Arguments
    ----------
    filename : string
        Filename with FASTQC data
    backend : string
        Database backend. Only this is required for an sqlite database.
    host : string
        Database host name
    username : string
        Database user name
    password : string
        Database password
    port : int
        Database server port.
    '''

    parser = CSV2DB.buildParser()
    (options, args) = parser.parse_args([])

    options.database_backend = backend
    options.database_host = host
    options.database_name = database
    options.database_username = username
    options.database_password = password
    options.database_port = port
    options.allow_empty = True

    for fn in glob.glob(filename):
        prefix = os.path.basename(os.path.dirname(fn))
        results = []

        for name, status, header, data in FastqcSectionIterator(
                IOTools.openFile(fn)):
            # do not collect basic stats, see loadFastQCSummary
            if name == "Basic Statistics":
                continue

            options.tablename = prefix + "_" + re.sub(" ", "_", name)

            inf = StringIO("\n".join([header] + data) + "\n")
            CSV2DB.run(inf, options)
            results.append((name, status))

        # load status table
        options.tablename = prefix + "_status"

        inf = StringIO(
            "\n".join(["name\tstatus"] +
                      ["\t".join(x) for x in results]) + "\n")
        CSV2DB.run(inf, options)


def buildFastQCSummaryStatus(infiles, outfile, datadir):
    '''collect fastqc status results from multiple runs into a single table.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    '''

    outf = IOTools.openFile(outfile, "w")
    names = set()
    results = []
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".fastqc")
        filename = os.path.join(datadir,
                                track + "*_fastqc",
                                "fastqc_data.txt")
        # there can be missing sections
        for fn in glob.glob(filename):
            stats = collections.defaultdict(str)
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                stats[name] = status

            results.append((track, fn, stats))
            names.update(list(stats.keys()))

    names = list(names)
    outf.write("track\tfilename\t%s\n" % "\t".join(names))
    for track, fn, stats in results:
        outf.write("%s\t%s\t%s\n" %
                   (track, os.path.dirname(fn),
                    "\t".join(stats[x] for x in names)))
    outf.close()


def buildFastQCSummaryBasicStatistics(infiles, outfile, datadir):
    '''collect fastqc summary results from multiple runs into a single table.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    '''

    data = collectFastQCSections(infiles, "Basic Statistics", datadir)

    outf = IOTools.openFile(outfile, "w")
    first = True
    for track, status, header, rows in data:
        rows = [x.split("\t") for x in rows]
        if first:
            headers = [row[0] for row in rows]
            outf.write("track\t%s\n" % "\t".join(headers))
            first = False
        outf.write("%s\t%s\n" % (track, "\t".join([row[1] for row in rows])))
    outf.close()


def buildExperimentReadQuality(infiles, outfile, datadir):
    """build per-experiment read quality summary.

    Arguments
    ---------
    infiles : list
        List of filenames with fastqc output (logging information). The
        track name is derived from that.
    outfile : list
        Output filename in :term:`tsv` format.
    datadir : string
        Location of actual Fastqc output to be parsed.

    """
    data = collectFastQCSections(infiles,
                                 "Per sequence quality scores",
                                 datadir)
    first = True

    if len(data) == 0:
        raise ValueError("received no data")

    for track, status, header, rows in data:
        T = track
        rows = [list(map(float, x.split("\t"))) for x in rows]
        header = header.split("\t")
        if first:
            first = False
            df_out = pd.DataFrame(rows)
            df_out.columns = header
            df_out.rename(columns={"Count": track}, inplace=True)
        else:
            df = pd.DataFrame(rows)
            df.columns = header
            df.rename(columns={"Count": track}, inplace=True)
            df_out = df_out.merge(df, how="outer", on="Quality", sort=True)

    df_out.set_index("Quality", inplace=True)
    df_out = pd.DataFrame(df_out.sum(axis=1))
    df_out.columns = ["_".join(T.split("-")[:-1]), ]

    df_out.to_csv(IOTools.openFile(outfile, "w"), sep="\t")
