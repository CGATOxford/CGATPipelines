'''cgat_logfiles2tsv.py - create summary from logfiles
===================================================


Purpose
-------

This script takes a list of logfiles and collates summary information
about execution times. This can be useful for post-mortem
benchmark analysis.

This script uses the ``# job finished`` tag that is added by scripts
using the module :mod:`CGAT.Experiment`.

Usage
-----

To collect logfile information from all files matching the pattern
``bwa.dir/mC-juvenile-stressed-R[12]*.log``, type::

   python cgat_logfiles2tsv.py --glob="bwa.dir/mC-juvenile-stressed-R[12]*.log"

to receive output such as this::

  file    chunks  wall    user    sys     cuser   csys
  bwa.dir/mC-juvenile-stressed-R2.bwa.bam.log     2       2552.00 1563.91  13.73    0.00    0.04
  bwa.dir/mC-juvenile-stressed-R2.bwa.bw.log      1       2068.00 170.66    4.50  237.51  1194.92
  bwa.dir/mC-juvenile-stressed-R1.bwa.bam.log     2       1378.00 762.52    9.90    0.00    0.04
  bwa.dir/mC-juvenile-stressed-R1.bwa.contextstats.log    1       948.00  150.21    2.13  726.00    7.92
  bwa.dir/mC-juvenile-stressed-R2.bwa.contextstats.log    1       935.00  137.00    2.26  775.07    8.35
  bwa.dir/mC-juvenile-stressed-R1.bwa.bw.log      1       2150.00 159.64    4.12  214.59  1566.41
  total   8       10031.00        2943.94  36.64  1953.17 2777.68

The output lists for each file how often it was executed (``chunks``) and
the total execution time in terms of wall clock time, user time, system
time, child process user time and child process system time.

The last line contains the sum total.

Type::

   python cgat_logfiles2tsv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import gzip
import glob

import CGAT.Experiment as E
import CGAT.Logfile as Logfile


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-g", "--glob", dest="glob_pattern", type="string",
        help="glob pattern to use for collecting files [%default].")

    parser.add_option(
        "-f", "--file-pattern", dest="file_pattern", type="string",
        help="only check files matching this pattern [%default].")

    parser.add_option("-m", "--mode", dest="mode", type="choice",
                      choices=("file", "node"),
                      help="analysis mode [%default].")

    parser.add_option(
        "-r", "--recursive", action="store_true",
        help="recursively look for logfiles from current directory "
        "[%default].")

    parser.set_defaults(
        truncate_sites_list=0,
        glob_pattern="*.log",
        mode="file",
        recursive=False,
    )

    (options, args) = E.Start(parser)

    if args:
        filenames = args
    elif options.glob_pattern:
        filenames = glob.glob(options.glob_pattern)

    if len(filenames) == 0:
        raise ValueError("no files to analyse")

    if options.mode == "file":
        totals = Logfile.LogFileData()

        options.stdout.write("file\t%s\n" % totals.getHeader())

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open(filename, "r")
            else:
                infile = open(filename, "r")

            subtotals = Logfile.LogFileData()
            for line in infile:
                subtotals.add(line)

            infile.close()

            options.stdout.write("%s\t%s\n" % (filename, str(subtotals)))
            totals += subtotals

        options.stdout.write("%s\t%s\n" % ("total", str(totals)))

    elif options.mode == "node":

        chunks_per_node = {}

        rx_node = re.compile("# job started at .* \d+ on (\S+)")

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open(filename, "r")
            else:
                infile = open(filename, "r")

            data = Logfile.LogFileDataLines()

            for line in infile:

                if rx_node.match(line):
                    node_id = rx_node.match(line).groups()[0]
                    data = Logfile.LogFileDataLines()
                    if node_id not in chunks_per_node:
                        chunks_per_node[node_id] = []
                    chunks_per_node[node_id].append(data)
                    continue

                data.add(line)

        options.stdout.write("node\t%s\n" % data.getHeader())
        total = Logfile.LogFileDataLines()

        for node, data in sorted(chunks_per_node.items()):
            subtotal = Logfile.LogFileDataLines()
            for d in data:
                # options.stdout.write( "%s\t%s\n" % (node, str(d) ) )
                subtotal += d

            options.stdout.write("%s\t%s\n" % (node, str(subtotal)))

            total += subtotal

        options.stdout.write("%s\t%s\n" % ("total", str(total)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
