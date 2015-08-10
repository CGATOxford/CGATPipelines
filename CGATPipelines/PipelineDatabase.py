'''
PipelineDatabase.py - utility functions for working with a database
===================================================================

'''

import os
import CGATPipelines.Pipeline as P

# set from calling module
PARAMS = {}


def importFromIterator(
        outfile,
        tablename,
        iterator,
        columns=None,
        indices=None):
    '''import data in *iterator* into *tablename* via temporary file.

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
