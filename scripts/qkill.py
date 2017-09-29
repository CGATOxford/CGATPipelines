"""qkill.py - kill jobs in the queue
==================================


Purpose
-------

This script kills jobs in the job queue according to certain criteria.
It is complementary to the ``qdel`` command, but provides more
flexible pattern matching.

This script requires the ``qstat`` and ``qdel`` commands to be in the
user's PATH. The script works with Sun Grid Engine.

Usage
-----

For example, to kill all jobs matching the pattern ``pipeline``,
type::

  python qkill.py -p=pipeline

Type::

   python qkill.py --help

for command line help.

Command line options
--------------------

"""

import sys
import re
import subprocess
import CGAT.Experiment as E
import xml.etree.ElementTree
import io as StringIO


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-p", "--pattern-identifier", dest="pattern", type="string",
        help="jobs matching `pattern` in their job "
        "description will be killed [default=%default].")

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="do dry run, do not kill [default=%default].")

    parser.set_defaults(
        pattern=None,
        dry_run=False,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    output = StringIO.StringIO(
        subprocess.Popen(["qstat", "-xml"],
                         stdout=subprocess.PIPE).communicate()[0])

    tree = xml.etree.ElementTree.ElementTree(file=output)

    ntested = 0
    to_kill = set()

    if options.pattern:
        pattern = re.compile(options.pattern)
    else:
        pattern = None

    for x in tree.getiterator("job_list"):
        ntested += 1
        id = x.find("JB_job_number").text
        name = x.find("JB_name").text
        if pattern and pattern.search(name):
            to_kill.add(id)

    nkilled = len(to_kill)
    if not options.dry_run:
        p = subprocess.Popen(
            ["qdel", ",".join(to_kill)], stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()

    E.info("ntested=%i, nkilled=%i" % (ntested, nkilled))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
