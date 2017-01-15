'''cgat_cwd2list.py - list directory contents
========================================

:Date: |today|
:Tags: Python

Purpose
-------

create a flat file that lists the current contents of the current
directory and its subfolders.  Useful for when data is deleted due to
space contraints. If files need to be recreated then the difference
between the current state and previous states can be assessed.

Usage
-----

The command::

   python cgat_cwd2list.py

This will create a file such as :file:`CWD_2015-09-15_12:45:29`
listing all the files in the current directory and its subdirectories::

  ##contents of cwd on 2015-09-15_12:45:29

  ./scripts/cgat_cluster_distribute.rst
  ./scripts/cgat_tsv2links.rst
  ./scripts/qkill.rst
  ./scripts/farm.rst
  ./scripts/cgat_ruffus_profile.rst
  ./scripts/__init__.rst
  ./scripts/submit.rst
  ./scripts/cgat_clean.rst
  ...

Type::

   python cgat_cwd2list.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import time
import datetime

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    dir2files = {}
    for root, directory, files in os.walk("."):
        dir2files[root] = files

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')
    filename = "CWD_%s" % st
    E.info("outputting directory state to %s" % filename)
    with IOTools.openFile(filename, "w") as outf:
        outf.write("##contents of cwd on %s\n\n" % st)
        for directory, files in dir2files.items():
            for file in files:
                path = os.path.join(directory, file)
                outf.write(path + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
