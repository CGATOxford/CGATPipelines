'''nofarm.py - stand-in for ``farm.py``
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script is a stand-in for :doc:`farm`. Instead of splitting the
input and running jobs on the cluster, the command is executed as is.

This script accepts all the options of ``farm.py``, but ignores most
of them. This is useful to turn of remote-processing when debugging a
pipeline locally.

Usage
-----

Example::

   python nofarm.py --help

Type::

   python nofarm.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import subprocess
import farm

import CGAT.Experiment as E


def main(argv=None):

    parser = farm.getOptionParser()

    (options, args) = E.Start(parser,
                              add_cluster_options=True)

    cmd = args[0]
    if len(args) > 1:
        cmd += " '" + "' '".join(args[1:]) + "'"

    cmd = re.sub("%DIR%", "", cmd)
    retcode = subprocess.call(cmd,
                              shell=True,
                              stdin=sys.stdin,
                              stdout=sys.stdout,
                              cwd=os.getcwd(),
                              close_fds=True)
    E.Stop()

if __name__ == '__main__':
    main()
