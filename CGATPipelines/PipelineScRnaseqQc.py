import os
import subprocess
import CGAT.Experiment as E
import sqlite3 as sql
import pandas as pd
import pandas.io.sql as pdsql
import re
import numpy as np

# ------------------------------------------------------- #
# Functions for expression quantification
# ------------------------------------------------------- #


def runSailfishIndex(fasta_file, outdir, threads,
                     kmer):
    '''
    Wrapper for sailfish index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fasta_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = '''
    sailfish index --transcripts %s --out %s --threads %i --kmerSize %i
    ''' % (fasta_file, outdir, threads, kmer)

    os.system(command)


def runSailfishQuant(fasta_index, fastq_files, output_dir,
                     paired=False, library="ISF", threads=4,
                     gene_gtf=None):
    '''
    Wrapper for sailfish quant command
    '''

    decompress = False
    if len(fastq_files) > 1:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass
    else:
        if fastq_files[0].endswith(".gz"):
            decompress = True
        else:
            pass

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " sailfish quant --index %s -l %s  -o %s " % (fasta_index,
                                                            library,
                                                            output_dir)

    states.append(command)

    if threads:
        states.append(" --threads %i " % threads)
    else:
        pass

    if gene_gtf:
        states.append(" --geneMap %s " % gene_gtf)
    else:
        pass

    # sailfish does not handle compress files natively,
    # need to decompress on the fly with advanced
    # bash syntax
    if decompress and paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1.gz",
                                                                   fq)])
        fstr_format = " ".join(["%s" for hq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2.gz",
                                                                    sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 <(zcat %s)" % sdecomp_format

        states.append(decomp_second)

    elif decompress and not paired:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1.gz",
                                                                   fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -r <(zcat %s)" % fdecomp_format

        states.append(decomp_first)

    elif paired and not decompress:
        first_mates = tuple([fq for fq in fastq_files if re.search("fastq.1",
                                                                   fq)])
        fstr_format = " ".join(["%s" for sq in first_mates])
        fdecomp_format = fstr_format % first_mates
        decomp_first = " -1 %s " % fdecomp_format

        states.append(decomp_first)

        second_mates = tuple([sq for sq in fastq_files if re.search("fastq.2",
                                                                    sq)])
        sstr_format = " ".join(["%s" for aq in second_mates])
        sdecomp_format = sstr_format % second_mates
        decomp_second = " -2 %s " % sdecomp_format

        states.append(decomp_second)

    statement = " ".join(states)

    # subprocess cannot handle process substitution
    # therefore needs to be wrapped in /bin/bash -c '...'
    # for bash to interpret the substitution correctly

    process = subprocess.Popen(statement, shell=True,
                               executable="/bin/bash")

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "-------------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was \n%s\n%s\n"
            "-------------------------------------------" %
            (-process.returncode, stderr, statement))


def runKallistoIndex(fasta_file, outfile, kmer=31):
    '''
    Wrapper for kallisto index
    '''

    if fasta_file.endswith(".fa"):
        pass
    elif fast_file.endswith(".fasta"):
        pass
    else:
        E.warn("are you sure this is a fasta file?")

    command = "kallisto index --index=%s  %s" % (outfile,
                                                 fasta_file)

    os.system(command)


def runKallistoQuant(fasta_index, fastq_files, output_dir,
                     bias=False, bootstrap=None,
                     seed=1245, threads=None, plaintext=False):
    '''
    Wrapper for kallisto quant command
    '''

    if len(fastq_files) > 1:
        fastqs = " ".join(fastq_files)
    else:
        fastqs = fastq_files

    # check output directory is an absolute path
    if os.path.isabs(output_dir):
        pass
    else:
        out_dir = os.path.abspath(output_dir)

    states = []
    command = " kallisto quant --index=%s --output-dir=%s" % (fasta_index,
                                                              output_dir)
    states.append(command)

    if bias:
        states.append(" --use-bias ")
    else:
        pass

    if bootstrap:
        states.append(" --bootstrap=%i --seed=%i " % (bootstrap,
                                                      seed))
    else:
        pass

    if plaintext:
        states.append(" --plaintext ")
    else:
        pass

    if threads:
        states.append(" --threads=%i " % threads)
    else:
        pass

    states.append(" %s " % fastqs)

    statement = " ".join(states)

    # need to rename output files to conform to input/output
    # pattern as required.  Default name is abundance*.txt
    # when using plaintext output
    # kaliisto requires an output directory - create many small
    # directories, one for each file.
    # then extract the abundance.txt file and rename using the
    # input/output pattern

    os.system(statement)

# ----------------------------------------------------------- #
# miscellaneous/utility functions
# ----------------------------------------------------------- #


def getTableFromDb(db, table):
    '''
    Get a table from a database with pandas
    '''

    state = ''' SELECT * FROM %(table)s;''' % locals()
    dbh = sql.connect(db)
    df = pdsql.read_sql(state, dbh)
    df.index = df["track"]
    df.drop(labels="track", inplace=True,
            axis=1)

    return df


def cleanStatsTable(stats_file):
    '''
    Take in a table containing aggregated stats
    and clean by removing duplicate columns
    '''

    _df = pd.read_table(stats_file, sep="\t", header=0,
                        index_col=None, mangle_dupe_cols=False)
    # drop duplicates is case sensitive, convert all to
    # same case - SQL is not case sensitive so will throw
    # a hissy fit for same column names in different cases
    _df.columns = [cx.lower() for cx in _df.columns]
    _df = _df.T.drop_duplicates().T
    _df.index = _df["track"]
    return _df


def extractTranscriptCounts(con, table):
    '''
    Extract transcript model counts for a
    given sample

    Arguments
    ---------
    con: sqlite.connection
      An SQLite connection

    table: string
      the table to extract the transcript counts
      from.

    Returns
    -------
    coverages: pandas.Core.Series
    '''

    statement = '''
    SELECT coverage_sense_pcovered
    FROM %(table)s
    WHERE coverage_sense_nval > 0;
    ''' % locals()

    coverages = pdsql.read_sql(statement, con)
    coverages = coverages.loc[:, "coverage_sense_pcovered"]
    return coverages


def summariseOverBins(coverages, bins):
    '''
    Summarise model coverages over a set of bins

    Argumnets
    ---------
    coverages: pandas.Core.Series
      coverages over gene/transcripts

    bins: list
      values corresponding to percentage bins

    Returns
    -------
    freqs: numpy.array
      frequency array of coverages over percentiles
    '''

    freqs = np.zeros(shape=len(bins),
                     dtype=np.float64)
    for i in range(len(bins)):
        if i == 0:
            hits = coverages <= bins[i]
        else:
            hits = (coverages <= bins[i]) & (coverages > bins[i - 1])

        freqs[i] = len(coverages[hits])

    return freqs


def getModelCoverage(db, table_regex, model_type="transcript"):
    '''
    Compute transcript model coverage stats

    Arguments
    ---------
    db: string
      database containing transcript counts

    table_regex: string
      regular expression for transcript count table

    model_type: string
      calculate coverages over either transcripts or
      genes.  Default is gene models

    Returns
    -------
    coverage_df: Pandas.Core.DataFrame
      model coverage stats summarised for each cell
    '''

    # need to regex for all the tables, one for each sample
    # fetch_all returns a list of tuples
    dbh = sql.connect(db)
    cursor = dbh.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

    tab_reg = re.compile(table_regex)
    table_list = [tx[0] for tx in cursor.fetchall() if re.search(tab_reg,
                                                                 tx[0])]

    # pull out counts for each cell and compute coverages

    bins = list(range(0, 101))
    cov_dict = {}
    for tab in table_list:
        covs = extractTranscriptCounts(dbh, tab)
        freq_array = summariseOverBins(covs, bins)
        cov_dict[tab] = freq_array

    coverage_df = pd.DataFrame(cov_dict).T
    # create a regex group to remove superfluous characters
    # from the track names
    ix_re = re.compile(
        "_(?P<run>\d+)_(?P<plate>\d+)_(?P<well>\d+)_(?P<mapper>\S+)_transcript_counts")
    re_matches = [re.match(ix_re, ix) for ix in coverage_df.index]
    indx = ["%s_%s-%s.%s" % rm.group(1, 2, 3, 4) for rm in re_matches]
    coverage_df.index = indx
    coverage_df.columns = ["Bin%i" % bx for bx in coverage_df.columns]
    return coverage_df
