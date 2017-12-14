"""Database.py - Database upload for ruffus pipelines
=========================================================

Reference
---------

"""
import re
import os
import sqlite3
from CGAT import Database as Database
import CGAT.Experiment as E

from CGAT.IOTools import touchFile, snip

from CGATPipelines.Pipeline.Execution import buildStatement, run
from CGATPipelines.Pipeline.Files import getTempFile
from CGATPipelines.Pipeline.Parameters import getParams


def tablequote(track):
    '''quote a track name such that is suitable as a table name.'''
    return re.sub("[-(),\[\].]", "_", track)


def toTable(outfile):
    '''convert a filename from a load statement into a table name.

    This method checks if the filename ends with ".load". The suffix
    is then removed and the filename quoted so that it is suitable
    as a table name.

    Arguments
    ---------
    outfile : string
        A filename ending in ".load".

    Returns
    -------
    tablename : string

    '''
    assert outfile.endswith(".load")
    name = os.path.basename(outfile[:-len(".load")])
    return tablequote(name)


def build_load_statement(tablename, retry=True, options=""):
    """build a command line statement to upload data.

    Upload is performed via the :doc:`csv2db` script.

    The returned statement is suitable to use in pipe expression.
    This method is aware of the configuration values for database
    access and the chosen database backend.

    For example::

        load_statement = P.build_load_statement("data")
        statement = "cat data.txt | %(load_statement)s"
        P.run()

    Arguments
    ---------
    tablename : string
        Tablename for upload
    retry : bool
        Add the ``--retry`` option to `csv2db.py`
    options : string
        Command line options to be passed on to `csv2db.py`

    Returns
    -------
    string

    """

    opts = []

    if retry:
        opts.append(" --retry ")

    PARAMS = getParams()
    backend = PARAMS["database_backend"]

    if backend not in ("sqlite", "mysql", "postgres"):
        raise NotImplementedError(
            "backend %s not implemented" % backend)

    opts.append("--database-backend=%s" % backend)
    opts.append("--database-name=%s" %
                PARAMS.get("database_name"))
    opts.append("--database-host=%s" %
                PARAMS.get("database_host", ""))
    opts.append("--database-user=%s" %
                PARAMS.get("database_username", ""))
    opts.append("--database-password=%s" %
                PARAMS.get("database_password", ""))
    opts.append("--database-port=%s" %
                PARAMS.get("database_port", 3306))

    db_options = " ".join(opts)

    statement = ('''
    cgat csv2db
    %(db_options)s
    %(options)s
    --table=%(tablename)s
    ''')

    load_statement = buildStatement(**locals())

    return load_statement


def load(infile,
         outfile=None,
         options="",
         collapse=False,
         transpose=False,
         tablename=None,
         retry=True,
         limit=0,
         shuffle=False,
         job_memory=None):
    """import data from a tab-separated file into database.

    The table name is given by outfile without the
    ".load" suffix.

    A typical load task in ruffus would look like this::

        @transform("*.tsv.gz", suffix(".tsv.gz"), ".load")
        def loadData(infile, outfile):
            P.load(infile, outfile)

    Upload is performed via the :doc:`csv2db` script.

    Arguments
    ---------
    infile : string
        Filename of the input data
    outfile : string
        Output filename. This will contain the logging information. The
        table name is derived from `outfile` if `tablename` is not set.
    options : string
        Command line options for the `csv2db.py` script.
    collapse : string
        If set, the table will be collapsed before loading. This
        transforms a data set with two columns where the first column
        is the row name into a multi-column table.  The value of
        collapse is the value used for missing values.
    transpose : string
        If set, the table will be transposed before loading. The first
        column in the first row will be set to the string within
        transpose.
    retry : bool
        If True, multiple attempts will be made if the data can
        not be loaded at the first try, for example if a table is locked.
    limit : int
        If set, only load the first n lines.
    shuffle : bool
        If set, randomize lines before loading. Together with `limit`
        this permits loading a sample of rows.
    job_memory : string
        Amount of memory to allocate for job. If unset, uses the global
        default.
    """
    PARAMS = getParams()
    if job_memory is None:
        job_memory = PARAMS["cluster_memory_default"]

    if not tablename:
        tablename = toTable(outfile)

    statement = []

    if infile.endswith(".gz"):
        statement.append("zcat %(infile)s")
    else:
        statement.append("cat %(infile)s")

    if collapse:
        statement.append(
            "cgat table2table --collapse=%(collapse)s")

    if transpose:
        statement.append(
            """cgat table2table --transpose
            --set-transpose-field=%(transpose)s""")

    if shuffle:
        statement.append("cgat randomize_lines --keep-header=1")

    if limit > 0:
        # use awk to filter in order to avoid a pipeline broken error from head
        statement.append("awk 'NR > %i {exit(0)} {print}'" % (limit + 1))
        # ignore errors from cat or zcat due to broken pipe
        ignore_pipe_errors = True

    statement.append(build_load_statement(tablename,
                                          options=options,
                                          retry=retry))

    statement = " | ".join(statement) + " > %(outfile)s"

    to_cluster = False

    run()


def concatenateAndLoad(infiles,
                       outfile,
                       regex_filename=None,
                       header=None,
                       cat="track",
                       has_titles=True,
                       missing_value="na",
                       retry=True,
                       tablename=None,
                       options="",
                       job_memory=None):
    """concatenate multiple tab-separated files and upload into database.

    The table name is given by outfile without the
    ".load" suffix.

    A typical concatenate and load task in ruffus would look like this::

        @merge("*.tsv.gz", ".load")
        def loadData(infile, outfile):
            P.concatenateAndLoad(infiles, outfile)

    Upload is performed via the :doc:`csv2db` script.

    Arguments
    ---------
    infiles : list
        Filenames of the input data
    outfile : string
        Output filename. This will contain the logging information. The
        table name is derived from `outfile`.
    regex_filename : string
        If given, *regex_filename* is applied to the filename to extract
        the track name. If the pattern contains multiple groups, they are
        added as additional columns. For example, if `cat` is set to
        ``track,method`` and `regex_filename` is ``(.*)_(.*).tsv.gz``
        it will add the columns ``track`` and method to the table.
    header : string
        Comma-separated list of values for header.
    cat : string
        Column title for column containing the track name. The track name
        is derived from the filename, see `regex_filename`.
    has_titles : bool
        If True, files are expected to have column titles in their first row.
    missing_value : string
        String to use for missing values.
    retry : bool
        If True, multiple attempts will be made if the data can
        not be loaded at the first try, for example if a table is locked.
    tablename: string
        Name to use for table. If unset derive from outfile.
    options : string
        Command line options for the `csv2db.py` script.
    job_memory : string
        Amount of memory to allocate for job. If unset, uses the global
        default.

    """
    PARAMS = getParams()

    if job_memory is None:
        job_memory = PARAMS["cluster_memory_default"]

    if tablename is None:
        tablename = toTable(outfile)

    infiles = " ".join(infiles)

    passed_options = options
    load_options, cat_options = ["--add-index=track"], []

    if regex_filename:
        cat_options.append("--regex-filename='%s'" % regex_filename)

    if header:
        load_options.append("--header-names=%s" % header)

    if not has_titles:
        cat_options.append("--no-titles")

    cat_options = " ".join(cat_options)
    load_options = " ".join(load_options) + " " + passed_options

    load_statement = build_load_statement(tablename,
                                          options=load_options,
                                          retry=retry)

    statement = '''cgat combine_tables
    --cat=%(cat)s
    --missing-value=%(missing_value)s
    %(cat_options)s
    %(infiles)s
    | %(load_statement)s
    > %(outfile)s'''

    to_cluster = False

    run()


def mergeAndLoad(infiles,
                 outfile,
                 suffix=None,
                 columns=(0, 1),
                 regex=None,
                 row_wise=True,
                 retry=True,
                 options="",
                 prefixes=None):
    '''merge multiple categorical tables and load into a database.

    The tables are merged and entered row-wise, i.e, the contents of
    each file are a row.

    For example, the statement::

        mergeAndLoad(['file1.txt', 'file2.txt'],
                     "test_table.load")

    with the two files::
        > cat file1.txt
        Category    Result
        length      12
        width       100

        > cat file2.txt
        Category    Result
        length      20
        width       50

    will be added into table ``test_table`` as::
        track   length   width
        file1   12       100
        file2   20       50

    If row-wise is set::
        mergeAndLoad(['file1.txt', 'file2.txt'],
                     "test_table.load", row_wise=True)

    ``test_table`` will be transposed and look like this::
        track    file1 file2
        length   12    20
        width    20    50

    Arguments
    ---------
    infiles : list
        Filenames of the input data
    outfile : string
        Output filename. This will contain the logging information. The
        table name is derived from `outfile`.
    suffix : string
        If `suffix` is given, the suffix will be removed from the filenames.
    columns : list
        The columns to be taken. By default, the first two columns are
        taken with the first being the key. Filenames are stored in a
        ``track`` column. Directory names are chopped off.  If
        `columns` is set to None, all columns will be taken. Here,
        column names will receive a prefix given by `prefixes`. If
        `prefixes` is None, the filename will be added as a prefix.
    regex : string
        If set, the full filename will be used to extract a
        track name via the supplied regular expression.
    row_wise : bool
        If set to False, each table will be a column in the resulting
        table.  This is useful if histograms are being merged.
    retry : bool
        If True, multiple attempts will be made if the data can
        not be loaded at the first try, for example if a table is locked.
    options : string
        Command line options for the `csv2db.py` script.
    prefixes : list
        If given, the respective prefix will be added to each
        column. The number of `prefixes` and `infiles` needs to be the
        same.

    '''
    PARAMS = getParams()
    if len(infiles) == 0:
        raise ValueError("no files for merging")

    if suffix:
        header = ",".join([os.path.basename(snip(x, suffix)) for x in infiles])
    elif regex:
        header = ",".join(["-".join(re.search(regex, x).groups())
                          for x in infiles])
    else:
        header = ",".join([os.path.basename(x) for x in infiles])

    header_stmt = "--header-names=%s" % header

    if columns:
        column_filter = "| cut -f %s" % ",".join(map(str,
                                                 [x + 1 for x in columns]))
    else:
        column_filter = ""
        if prefixes:
            assert len(prefixes) == len(infiles)
            header_stmt = "--prefixes=%s" % ",".join(prefixes)
        else:
            header_stmt = "--add-file-prefix"

    if infiles[0].endswith(".gz"):
        filenames = " ".join(
            ["<( zcat %s %s )" % (x, column_filter) for x in infiles])
    else:
        filenames = " ".join(
            ["<( cat %s %s )" % (x, column_filter) for x in infiles])

    if row_wise:
        transform = """| perl -p -e "s/bin/track/"
        | cgat table2table --transpose""" % PARAMS
    else:
        transform = ""

    load_statement = build_load_statement(
        toTable(outfile),
        options="--add-index=track " + options,
        retry=retry)

    statement = """cgat combine_tables
    %(header_stmt)s
    --skip-titles
    --missing-value=0
    --ignore-empty
    %(filenames)s
    %(transform)s
    | %(load_statement)s
    > %(outfile)s
    """

    to_cluster = False

    run()


def connect():
    """connect to SQLite database used in this pipeline.

    .. note::
       This method is currently only implemented for sqlite
       databases. It needs refactoring for generic access.
       Alternatively, use an full or partial ORM.

    If ``annotations_database`` is in PARAMS, this method
    will attach the named database as ``annotations``.

    Returns
    -------
    dbh
       a database handle

    """

    # Note that in the future this might return an sqlalchemy or
    # db.py handle.
    PARAMS = getParams()
    if PARAMS["database_backend"] == "sqlite":
        dbh = sqlite3.connect(getDatabaseName())

        if "annotations_database" in PARAMS:
            statement = '''ATTACH DATABASE '%s' as annotations''' % \
                        (PARAMS["annotations_database"])
            cc = dbh.cursor()
            cc.execute(statement)
            cc.close()
    else:
        raise NotImplementedError(
            "backend %s not implemented" % PARAMS["database_backend"])
    return dbh


def createView(dbhandle, tables, tablename, outfile,
               view_type="TABLE",
               ignore_duplicates=True):
    '''create a database view for a list of tables.

    This method performs a join across multiple tables and stores the
    result either as a view or a table in the database.

    Arguments
    ---------
    dbhandle :
        A database handle.
    tables : list of tuples
        Tables to merge. Each tuple contains the name of a table and
        the field to join with the first table. For example::

            tables = (
                "reads_summary", "track",
                "bam_stats", "track",
                "context_stats", "track",
                "picard_stats_alignment_summary_metrics", "track")

    tablename : string
        Name of the view or table to be created.
    outfile : string
        Output filename for status information.
    view_type : string
        Type of view, either ``VIEW`` or ``TABLE``.  If a view is to be
        created across multiple databases, use ``TABLE``.
    ignore_duplicates : bool
        If set to False, duplicate column names will be added with the
        tablename as prefix. The default is to ignore.

    '''

    Database.executewait(
        dbhandle,
        "DROP %(view_type)s IF EXISTS %(tablename)s" % locals())

    tracks, columns = [], []
    tablenames = [x[0] for x in tables]
    for table, track in tables:
        d = Database.executewait(
            dbhandle,
            "SELECT COUNT(DISTINCT %s) FROM %s" % (track, table))
        tracks.append(d.fetchone()[0])
        columns.append(
            [x.lower() for x in Database.getColumnNames(dbhandle, table)
             if x != track])

    E.info("creating %s from the following tables: %s" %
           (tablename, str(list(zip(tablenames, tracks)))))
    if min(tracks) != max(tracks):
        raise ValueError(
            "number of rows not identical - will not create view")

    from_statement = " , ".join(
        ["%s as t%i" % (y[0], x) for x, y in enumerate(tables)])
    f = tables[0][1]
    where_statement = " AND ".join(
        ["t0.%s = t%i.%s" % (f, x + 1, y[1])
         for x, y in enumerate(tables[1:])])

    all_columns, taken = [], set()
    for x, c in enumerate(columns):
        i = set(taken).intersection(set(c))
        if i:
            E.warn("duplicate column names: %s " % i)
            if not ignore_duplicates:
                table = tables[x][0]
                all_columns.extend(
                    ["t%i.%s AS %s_%s" % (x, y, table, y) for y in i])
                c = [y for y in c if y not in i]

        all_columns.extend(["t%i.%s" % (x, y) for y in c])
        taken.update(set(c))

    all_columns = ",".join(all_columns)
    statement = '''
    CREATE %(view_type)s %(tablename)s AS SELECT t0.track, %(all_columns)s
    FROM %(from_statement)s
    WHERE %(where_statement)s
    ''' % locals()

    Database.executewait(dbhandle, statement)

    nrows = Database.executewait(
        dbhandle, "SELECT COUNT(*) FROM view_mapping").fetchone()[0]

    if nrows == 0:
        raise ValueError(
            "empty view mapping, check statement = %s" %
            (statement % locals()))
    if nrows != min(tracks):
        E.warn("view creates duplicate rows, got %i, expected %i" %
               (nrows, min(tracks)))

    E.info("created view_mapping with %i rows" % nrows)
    touchFile(outfile)


def getDatabaseName():
    '''Return the database name associated with the pipeline.

    This method lookis in different sections in the ini file to permit
    both old style ``database`` and new style ``database_name``.

    This method has been implemented for backwards compatibility.

    Returns
    -------
    databasename : string
        Database name. Returns empty string if not found.

    Raises
    ------
    KeyError
       If no database name is found

    '''

    locations = ["database_name", "database"]
    PARAMS = getParams()
    for location in locations:
        database = PARAMS.get(location, None)
        if database is not None:
            return database

    raise KeyError("database name not found")


def importFromIterator(
        outfile,
        tablename,
        iterator,
        columns=None,
        indices=None):
    '''import data from an iterator into a database.

    Arguments
    ---------
    outfile : string
        Output file name
    tablename : string
        Table name
    iterator : iterator
        Iterator to import data from. The iterator should
        yield either list/tuples or dictionaries for each
        row in the table.
    columns : list
        Column names. If not given, the assumption is that
        iterator will dictionaries and column names are derived
        from that.
    indices : list
        List of column names to add indices on.
    '''

    tmpfile = getTempFile(".")

    if columns:
        keys, values = list(zip(*list(columns.items())))
        tmpfile.write("\t".join(values) + "\n")

    for row in iterator:
        if not columns:
            keys = list(row[0].keys())
            values = keys
            columns = keys
            tmpfile.write("\t".join(values) + "\n")

        tmpfile.write("\t".join(str(row[x]) for x in keys) + "\n")

    tmpfile.close()

    if indices:
        indices = " ".join("--add-index=%s" % x for x in indices)
    else:
        indices = ""

    load(tmpfile.name,
         outfile,
         tablename=tablename,
         options=indices)

    os.unlink(tmpfile.name)
