'''
cgat_check_deps.py - check whether the software dependencies are on your PATH
=============================================================================

Purpose
-------

.. This script takes the name of a CGAT pipeline and checks whether all
programs used in the pipeline (via "statement" and "P.run()") can be
located on the PATH.

Usage
-----

.. cgat cgat_check_deps --pipeline </path/to/pipeline_name.py>

Example::

   cgat cgat_check_deps --pipeline CGATPipelines/pipeline_annotations.py

Type::

   cgat cgat_check_deps --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import subprocess
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-p", "--pipeline", dest="pipeline", type="string",
                      help="Path to pipeline script")

    parser.add_option("-s", "--print-summary", dest="summary", action="store_true", default=False,
                      help="Print how many times a program is used [default=%default]")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # check existence of pipeline script
    if not os.access(options.pipeline, os.R_OK):
        raise IOError("Pipeline %s was not found\n" % options.pipeline)

    # parse pipeline script
    statement = '''awk '/statement =/,/P.run()/ {print}' %s
                   | grep -v '#'
                   | tr '\n' ' '
                   | tr '"' ' '
                   | sed 's/statement =//g'
                   | sed 's/P.run()//g'
                   '''
    # the code below has been adapted from the Execution module of Pipeline.py
    # cleaning up of statement
    # remove new lines and superfluous spaces and tabs
    statement = " ".join(re.sub("\t+", " ", statement).split("\n")).strip()
    if statement.endswith(";"):
        statement = statement[:-1]

    process = subprocess.Popen(statement % options.pipeline,
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    # get process output
    stdout, stderr = process.communicate()

    # check process status
    if process.returncode != 0:
        raise OSError(
            "Child was terminated by signal %i: \n"
            "The stderr was: \n%s\n%s\n" %
             (-process.returncode, stderr, statement))

    print stdout

    # dictionary where:
    # key = program name
    # value = number of times it has been called
    deps = {}

    for slice in stdout.split("'''"):
        # clean up white duplicated spaces
        clean_slice = ' '.join(slice.split())
        if len(clean_slice) > 0:
            for command in clean_slice.split("|"):
                # take program name
                groups  = re.match("^\s*(\w+)", command)
                if groups is not None:
                    # program name is first match
                    prog_name = groups.group(0)
                    # clean up duplicated white spaces
                    prog_name = ' '.join(prog_name.split())
                    if prog_name not in deps:
                        deps[prog_name] = 1
                    else:
                        deps[prog_name]+= 1

    check_path_failures = []

    # print dictionary ordered by value
    for k in sorted(deps, key=deps.get, reverse=True):
        if options.summary:
            print 'Program: {0!s} used {1} time(s) \n'.format(k,deps[k])
        if len(IOTools.which(k)) == 0:
            check_path_failures.append(k)
            print 'Program: {0!s} is NOT on your PATH!\n'.format(k)

    n_failures = len(check_path_failures)
    if n_failures == 0:
        print 'Congratulations! All required programs are available on your PATH\n'
    else:
        print 'The following programs are not on your PATH\n'
        for p in check_path_failure:
            print '{0!s}\n'.format(p)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

