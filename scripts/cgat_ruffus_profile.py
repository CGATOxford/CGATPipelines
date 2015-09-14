
'''cgat_ruffus_profile.py - analyze ruffus logfile
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script collects information about tasks that have completed or
are still running in a pipeline. It works by examining the logfile
:file:`pipeline.log` looking for the last active run. It will collect
a list of all tasks that have been executed or have just started and
display runtime information.

Usage
-----

To use the script, simply execute in a working directory of a pipeline::

  python cgat_ruffus_profile.py

Example output looks like this::

  section object  ncalls  duration        percall running
  task    runMedipsDMR    1       0        0.000  1
  task    buildWindowsFoldChangesPerMedian        1       4        4.000  0
  task    summarizeAllWindowsReadCounts   1       0        0.000  1
  task    loadPicardDuplicateStats        1       31      31.000  0
  task    (mkdir2)beforerunDESeq  2       11       5.500  0
  task    summarizeWindowsReadCounts      1       60274   60274.000       0
  task    runFilterAnalysis       1       0        0.000  1
  task    (mkdir1)beforerunEdgeR  2       7        3.500  0
  task    loadWindowsFoldChanges  1       126     126.000 0
  ...
  #//

  # running tasks
  runMedipsDMR
  summarizeAllWindowsReadCounts
  ...
  #//

  section object  ncalls  duration        percall running
  job     overlaps.dir/method_designJuvenilemCNoOutlier_mC-juvenile-stressed.overlap      1       37      37.000  0
  job     overlaps.dir/edger.overlap      1       0        0.000  1
  job     roi.dir/designJuvenileHippocampusmC.pvalue      2       209     104.500 0
  job     medips.dir/designJuvenileCortexmC.qvalue        1       256     256.000 0
  job     overlaps.dir/by_method_designFoetalCortexmC.overlap     1       35      35.000  0
  job     dmr_pvalues.load        1       12      12.000  0
  ...
  #//

  # running jobs
  overlaps.dir/edger.overlap
  overlaps.dir/deseq_mC-foetal-sal.overlap
  overlaps.dir/method_designJuvenileCortexmC_indows.all.overlap
  medips.dir/designFoetalCortexmC.mergedwindows.all.bed.gz
  ...
  #//

The output is divided in four sections displaying information for
tasks and jobs. A task is a collection of jobs.  For each of the task
or job summary, there are two sections. The first section lists all
the tasks and indicates:

ncalls
    The number of time this task was executed. This corresponds to the
    number of jobs started in this particular task
duration
    The execution time total for this task in wall clock time (seconds).
    This is the time from when the task started to when it was complete.
percall
    The average execution time per job for this task.
running
    Flag to indicate if the task is still running.

This section is followed by a list of all tasks that are still
running.

Next follow two sections describing equivalent information for jobs.

This script relies on the :file:`pipeline.log` being in a consistent
state. This might not be case if a pipeline has been executed several
times simultaneously. Use the option ``--ignore-errors`` to get
around this.

By default, the script will only look at the last execution of a
pipeline, use ``--no-reset`` to gather information over the full
logfile.

Type::

  python cgat_ruffus_profile.py --help

to get command line help.

Command line options
--------------------

'''

import sys
import os
import re
import datetime
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


class Counter(object):

    '''

    This class stores the calls per source as calls can be made by
    different sources in any order.
    '''

    def __init__(self):
        self._durations = collections.defaultdict(list)
        self._started = collections.defaultdict(int)
        self._calls = collections.defaultdict(int)

    def add(self, started, dt, source=None):

        if self._started[source] != 0 and not started:
            self._durations[source].append(dt - self._started[source])
            self._started[source] = 0
        elif self._started[source] == 0 and started:
            self._calls[source] += 1
            self._started[source] = dt
        else:
            raise ValueError(
                """inconsistent time points for %s, has_started=%s, is_started=%s.
                Possibly two pipelines have been running concurrently.
                """ % (source,
                       self._started[source], started))

    def reset(self, source=None):
        '''reset last event.'''
        self._started[source] = None
        self._calls[source] = 0

    def getDuration(self):
        if self._durations:
            x = None
            for source, durations in self._durations.iteritems():
                if x is None:
                    x = durations[0]
                for y in durations[1:]:
                    x += y
            return x
        else:
            return datetime.timedelta()

    def getCalls(self):
        return sum(self._calls.values())

    def getRunning(self):
        '''get numbers of tasks unfinished or still running.'''
        return len([x for x, y in self._started.iteritems() if y != 0])

    duration = property(getDuration)
    calls = property(getCalls)
    running = property(getRunning)


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-l", "--logfile", dest="logfile", type="string",
                      help="name of logfile [default=%default]")

    parser.add_option("-t", "--time", dest="time", type="choice",
                      choices=("seconds", "milliseconds"),
                      help="time to show [default=%default]")

    parser.add_option(
        "--no-reset", dest="reset", action="store_false",
        help="do not reset counters when a new pipeline run started "
        "The default is to reset so that only the counts from the latest "
        "pipeline execution are show "
        "[default=%default]")

    parser.add_option(
        "-f", "--filter-method", dest="filter", type="choice",
        choices=("unfinished", "running", "completed", "all"),
        help="apply filter to output [default=%default]")

    parser.add_option(
        "-i", "--ignore-errors", dest="ignore_errors", action="store_true",
        help="ignore errors [default=%default]")

    parser.set_defaults(sections=[],
                        logfile="pipeline.log",
                        filter="all",
                        reset=True,
                        time="seconds")

    (options, args) = E.Start(parser, argv)

    rx = re.compile("^[0-9]+")

    if options.sections:
        profile_sections = options.sections
    else:
        profile_sections = ("task", "job")

    counts = {}
    for section in profile_sections:
        counts[section] = collections.defaultdict(Counter)

    rootpath = os.path.abspath(".")

    infile = IOTools.openFile(options.logfile)

    for line in infile:
        if not rx.match(line):
            continue
        data = line[:-1].split()
        if len(data) < 5:
            continue
        date, time, level, source = data[:4]

        if re.search("output generated by", line):
            if options.reset:
                E.info("resetting counts at line=%s" % line[:-1])
                for section in profile_sections:
                    counts[section] = collections.defaultdict(Counter)
            continue

        if not re.match("task\.", source):
            continue

        dt = datetime.datetime.strptime(
            " ".join((date, time)), "%Y-%m-%d %H:%M:%S,%f")

        msg = "".join(data[4:])

        started_task, completed_task, started_job, completed_job = \
            (None, None, None, None)

        if re.search("task.log_at_level.\d+Task=(\S+)", msg):
            checked_task = re.search(
                "task.log_at_level.\d+Task=(\S+)", msg).groups()[0]
        elif re.search("Job=\[(\S+)->(\S+)\]Missingfile[s]*\[(\S+)\]", msg):
            started_infiles, started_job, missing = re.search(
                "Job=\[(\S+)->(\S+)\]Missingfile[s]*\[(\S+)\]", msg).groups()
        elif re.search("Job=\[(\S+)->(\S+)\]Missingfile[s]*", msg):
            started_infiles, started_job = re.search(
                "Job=\[(\S+)->(\S+)\]Missingfile[s]*", msg).groups()
        elif re.search("Taskentersqueue=(\S+)", msg):
            started_task = re.search("Taskentersqueue=(\S+)", msg).groups()[0]
        elif re.search("Job=\[(\S+)->(\S+)\]completed", msg):
            completed_infiles, completed_job = re.search(
                "Job=\[(\S+)->(\S+)\]completed", msg).groups()
        elif re.search("CompletedTask=(\S+)", msg):
            completed_task = re.search("CompletedTask=(\S+)", msg).groups()[0]
        elif re.search("UptodateTask=(\S+)", msg):
            completed_task = re.search("UptodateTask=(\S+)", msg).groups()[0]
        else:
            continue

        try:
            if started_task:
                counts["task"][started_task].add(True, dt, started_task)
            elif completed_task:
                counts["task"][completed_task].add(False, dt, completed_task)
            elif started_job:
                counts["job"][started_job].add(True, dt, started_job)
            elif completed_job:
                counts["job"][completed_job].add(False, dt, completed_job)
            else:
                raise ValueError("unknown action")
        except ValueError, msg:
            if not options.ignore_errors:
                raise

    if options.time == "milliseconds":
        f = lambda d: d.seconds + d.microseconds / 1000
    elif options.time == "seconds":
        f = lambda d: d.seconds + d.microseconds / 1000000

    for section in profile_sections:
        options.stdout.write("\t".join(
            ("section", "object", "ncalls",
             "duration", "percall", "running")) + "\n")

        running = []
        for objct, c in counts[section].iteritems():

            # apply filters
            if options.filter in ("unfinished", "running") and c.running == 0:
                continue

            d = f(c.duration)
            if c.calls > 0:
                percall = "%6.3f" % (d / float(c.calls))
            else:
                percall = "na"

            options.stdout.write("\t".join(
                (map(str,
                     (section, objct,
                      c.calls,
                      d,
                      percall,
                      c.running,
                      )))) + "\n")

            running.extend([x for x, y in c._started.iteritems() if y != 0])

        options.stdout.write("#//\n\n")

        if running:
            options.stdout.write("# running %ss\n" % section)
            options.stdout.write("\n".join(map(str, running)) + "\n")
            options.stdout.write("#//\n\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())
