'''
PipelineDatabase.py - utility functions for working with a database
===================================================================

This module contains utility functions for working with a database.

.. note::

   Some functionality of :mod:`Pipeline` could go in here.

'''

import os
import CGATPipelines.Pipeline as P


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

    tmpfile = P.getTempFile(".")

    if columns:
        keys, values = zip(*columns.items())
        tmpfile.write("\t".join(values) + "\n")

    for row in iterator:
        if not columns:
            keys = row[0].keys()
            values = keys
            columns = keys
            tmpfile.write("\t".join(values) + "\n")

        tmpfile.write("\t".join(str(row[x]) for x in keys) + "\n")

    tmpfile.close()

    if indices:
        indices = " ".join("--add-index=%s" % x for x in indices)
    else:
        indices = ""

    P.load(tmpfile.name,
           outfile,
           tablename=tablename,
           options=indices)

    os.unlink(tmpfile.name)
