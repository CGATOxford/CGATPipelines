'''
cgatflow.py - Computational Genomics Analysis Workflows
=======================================================

:Tags: Genomics

To use a specific workflow, type::

    cgatflow <workflow> [workflow options] [workflow arguments]

For this message and a list of available keywords type::

    cgatflow --help

To get help for a specify workflow, type::

    cgatflow <workflow> --help
'''

import os
import sys
import re
import glob
import imp
import collections
import CGATPipelines


def printListInColumns(l, ncolumns):
    '''output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % 3 != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):

    argv = sys.argv

    # paths to look for pipelines:
    path = os.path.join(os.path.abspath(os.path.dirname(CGATPipelines.__file__)))
    relpath = os.path.abspath("../src")

    paths = [path, relpath]
    
    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        pipelines = []
        for path in paths:
            pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of all available pipelines is:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))
        return

    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)
    
    (file, pathname, description) = imp.find_module(pipeline, paths)
    module = imp.load_module(pipeline, file, pathname, description)
    # remove 'cgatflow' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
