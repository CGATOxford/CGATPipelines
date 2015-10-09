##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
"""Control.py - Command line control for ruffus pipelines
=========================================================

The functions :func:`writeConfigFiles`, :func:`clean`,
:func:`clonePipeline` and :func:`peekParameters` provide the
functionality for particular pipeline commands.

:class:`MultiLineFormatter` improves the formatting
of long log messages, while
:class:`LoggingFilterRabbitMQ` intercepts ruffus log
messages and sends event information to a rabbitMQ message exchange
for task process monitoring.

Reference
---------

"""

import inspect
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from cStringIO import StringIO
from multiprocessing.pool import ThreadPool

# talking to mercurial
import hgapi

# talking to RabbitMQ
try:
    import pika
    HAS_PIKA = True
except ImportError:
    HAS_PIKA = False

# talking to a cluster
try:
    import drmaa
    HAS_DRMAA = True
except RuntimeError:
    HAS_DRMAA = False

from ruffus import pipeline_printout_graph, pipeline_printout, \
    pipeline_run, ruffus_exceptions, task


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from CGAT import Requirements as Requirements

from CGATPipelines.Pipeline.Utils import isTest, getCaller, getCallerLocals
from CGATPipelines.Pipeline.Execution import execute, startSession,\
    closeSession
from CGATPipelines.Pipeline.Local import getProjectName, getPipelineName

# Set from Pipeline.py
PARAMS = {}

# global options and arguments - set but currently not
# used as relevant sections are entered into the PARAMS
# dictionary. Could be deprecated and removed.
GLOBAL_OPTIONS, GLOBAL_ARGS = None, None


def writeConfigFiles(path):
    '''create default configuration files in `path`.
    '''

    for dest in ("pipeline.ini", "conf.py"):
        src = os.path.join(path, dest)
        if os.path.exists(dest):
            E.warn("file `%s` already exists - skipped" % dest)
            continue

        if not os.path.exists(src):
            raise ValueError("default config file `%s` not found" % src)

        shutil.copyfile(src, dest)
        E.info("created new configuration file `%s` " % dest)


def clonePipeline(srcdir, destdir=None):
    '''clone a pipeline.

    Cloning entails creating a mirror of the source pipeline.
    Generally, data files are mirrored by linking. Configuration
    files and the pipeline database will be copied.

    Without modification of any files, building the cloned pipeline in
    `destdir` should not re-run any commands. However, on deleting
    selected files, the pipeline should run from the appropriate
    point.  Newly created files will not affect the original pipeline.

    Cloning pipelines permits sharing partial results between
    pipelines, for example for parameter optimization.

    Arguments
    ---------
    scrdir : string
        Source directory
    destdir : string
        Destination directory. If None, use the current directory.

    '''

    if destdir is None:
        destdir = os.path.curdir

    E.info("cloning pipeline from %s to %s" % (srcdir, destdir))

    copy_files = ("conf.py", "pipeline.ini", "csvdb")
    ignore_prefix = (
        "report", "_cache", "export", "tmp", "ctmp",
        "_static", "_templates")

    def _ignore(p):
        for x in ignore_prefix:
            if p.startswith(x):
                return True
        return False

    for root, dirs, files in os.walk(srcdir):

        relpath = os.path.relpath(root, srcdir)
        if _ignore(relpath):
            continue

        for d in dirs:
            if _ignore(d):
                continue
            dest = os.path.join(os.path.join(destdir, relpath, d))
            os.mkdir(dest)
            # touch
            s = os.stat(os.path.join(root, d))
            os.utime(dest, (s.st_atime, s.st_mtime))

        for f in files:
            if _ignore(f):
                continue

            fn = os.path.join(root, f)
            dest_fn = os.path.join(destdir, relpath, f)
            if f in copy_files:
                shutil.copyfile(fn, dest_fn)
            else:
                # realpath resolves links - thus links will be linked to
                # the original target
                os.symlink(os.path.realpath(fn),
                           dest_fn)


def clean(files, logfile):
    '''clean up files given by glob expressions.

    Files are cleaned up by zapping, i.e. the files are set to size
    0. Links to files are replaced with place-holders.

    Information about the original file is written to `logfile`.

    Arguments
    ---------
    files : list
        List of glob expressions of files to clean up.
    logfile : string
        Filename of logfile.

    '''
    fields = ('st_atime', 'st_blksize', 'st_blocks',
              'st_ctime', 'st_dev', 'st_gid', 'st_ino',
              'st_mode', 'st_mtime', 'st_nlink',
              'st_rdev', 'st_size', 'st_uid')

    dry_run = PARAMS.get("dryrun", False)

    if not dry_run:
        if not os.path.exists(logfile):
            outfile = IOTools.openFile(logfile, "w")
            outfile.write("filename\tzapped\tlinkdest\t%s\n" %
                          "\t".join(fields))
        else:
            outfile = IOTools.openFile(logfile, "a")

    c = E.Counter()
    for fn in files:
        c.files += 1
        if not dry_run:
            stat, linkdest = IOTools.zapFile(fn)
            if stat is not None:
                c.zapped += 1
                if linkdest is not None:
                    c.links += 1
                outfile.write("%s\t%s\t%s\t%s\n" % (
                    fn,
                    time.asctime(time.localtime(time.time())),
                    linkdest,
                    "\t".join([str(getattr(stat, x)) for x in fields])))

    E.info("zapped: %s" % (c))
    outfile.close()

    return c


def peekParameters(workingdir,
                   pipeline,
                   on_error_raise=None,
                   prefix=None,
                   update_interface=False,
                   restrict_interface=False):
    '''peek configuration parameters from external pipeline.

    As the paramater dictionary is built at runtime, this method
    executes the pipeline in workingdir, dumping its configuration
    values and reading them into a dictionary.

    If either `pipeline` or `workingdir` are not found, an error is
    raised. This behaviour can be changed by setting `on_error_raise`
    to False. In that case, an empty dictionary is returned.

    Arguments
    ---------
    workingdir : string
       Working directory. This is the directory that the pipeline
       was executed in.
    pipeline : string
       Name of the pipeline script. The pipeline is assumed to live
       in the same directory as the current pipeline.
    on_error_raise : Bool
       If set to a boolean, an error will be raised (or not) if there
       is an error during parameter peeking, for example if
       `workingdir` can not be found. If `on_error_raise` is None, it
       will be set to the default, which is to raise an exception
       unless the calling script is imported or the option
       ``--is-test`` has been passed at the command line.
    prefix : string
       Add a prefix to all parameters. This is useful if the paramaters
       are added to the configuration dictionary of the calling pipeline.
    update_interface : bool
       If True, this method will prefix any options in the
       ``[interface]`` section with `workingdir`. This allows
       transparent access to files in the external pipeline.
    restrict_interface : bool
       If  True, only interface parameters will be imported.

    Returns
    -------
    config : dict
        Dictionary of configuration values.

    '''
    caller_locals = getCallerLocals()

    # check if we should raise errors
    if on_error_raise is None:
        on_error_raise = not isTest() and \
            "__name__" in caller_locals and \
            caller_locals["__name__"] == "__main__"

    # patch - if --help or -h in command line arguments,
    # do not peek as there might be no config file.
    if "--help" in sys.argv or "-h" in sys.argv:
        return {}

    # Attempt to locate directory with pipeline source code. This is a
    # patch as pipelines might be called within the repository
    # directory or from an installed location
    dirname = PARAMS["pipelinedir"]

    # called without a directory, use current directory
    if dirname == "":
        dirname = os.path.abspath(".")
    else:
        # if not exists, assume we want version located
        # in directory of calling script.
        if not os.path.exists(dirname):
            # directory is path of calling script
            dirname = os.path.dirname(caller_locals['__file__'])

    pipeline = os.path.join(dirname, pipeline)
    if not os.path.exists(pipeline):
        if on_error_raise:
            raise ValueError(
                "can't find pipeline at %s" % (pipeline))
        else:
            return {}

    if workingdir == "":
        workingdir = os.path.abspath(".")

    # patch for the "config" target - use default
    # pipeline directory if directory is not specified
    # working dir is set to "?!"
    if "config" in sys.argv or "check" in sys.argv or "clone" in sys.argv and workingdir == "?!":
        workingdir = os.path.join(PARAMS.get("pipelinedir"),
                                  IOTools.snip(pipeline, ".py"))

    if not os.path.exists(workingdir):
        if on_error_raise:
            raise ValueError(
                "can't find working dir %s" % workingdir)
        else:
            return {}

    statement = "python %s -f -v 0 dump" % pipeline
    process = subprocess.Popen(statement,
                               cwd=workingdir,
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    # process.stdin.close()
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            ("Child was terminated by signal %i: \n"
             "The stderr was: \n%s\n") %
            (-process.returncode, stderr))

    dump = None
    for line in stdout.split("\n"):
        if line.startswith("dump"):
            exec(line)

    # update interface
    if update_interface:
        for key, value in dump.items():
            if key.startswith("interface"):
                dump[key] = os.path.join(workingdir, value)

    # keep only interface if so required
    if restrict_interface:
        dump = dict([(k, v) for k, v in dump.iteritems()
                     if k.startswith("interface")])

    # prefix all parameters
    if prefix is not None:
        dump = dict([("%s%s" % (prefix, x), y) for x, y in dump.items()])

    return dump


class MultiLineFormatter(logging.Formatter):
    """add identation for multi-line entries.
    """

    def format(self, record):
        s = logging.Formatter.format(self, record)
        if record.message:
            header, footer = s.split(record.message)
            s = s.replace('\n', '\n' + ' ' * len(header))
        return s


class LoggingFilterRabbitMQ(logging.Filter):
    """pass event information to a rabbitMQ message queue.

    This is a log filter which detects messages from ruffus_ and sends
    them to a rabbitMQ message queue.

    A :term:`task` is a ruffus_ decorated function, which will execute
    one or more :term:`jobs`.

    Valid task/job status:

    update
       task/job needs updating
    completed
       task/job completed successfully
    failed
       task/job failed
    running
       task/job is running
    ignore
       ignore task/job (is up-to-date)

    Arguments
    ---------
    ruffus_text : string
        Log messages from ruffus.pipeline_printout. These are used
        to collect all tasks that will be executed during pipeline
        executation.
    project_name : string
        Name of the project
    pipeline_name : string
        Name of the pipeline
    host : string
        RabbitMQ host name
    exchange : string
        RabbitMQ exchange name

    """

    def __init__(self, ruffus_text,
                 project_name,
                 pipeline_name,
                 host="localhost",
                 exchange="ruffus_pipelines"):

        self.project_name = project_name
        self.pipeline_name = pipeline_name
        self.exchange = exchange

        # dictionary of jobs to run
        self.jobs = {}
        self.tasks = {}

        if not HAS_PIKA:
            self.connected = False
            return

        def split_by_job(text):
            text = "".join(text)
            job_message = ""
            # ignore first entry which is the docstring
            for line in text.split(" Job  = ")[1:]:
                try:
                    # long file names cause additional wrapping and
                    # additional white-space characters
                    job_name = re.search(
                        "\[.*-> ([^\]]+)\]", line).groups()
                except AttributeError:
                    raise AttributeError("could not parse '%s'" % line)
                job_status = "ignore"
                if "Job needs update" in line:
                    job_status = "update"

                yield job_name, job_status, job_message

        def split_by_task(text):
            block, task_name = [], None
            task_status = None
            for line in text.split("\n"):
                if line.startswith("Tasks which will be run"):
                    task_status = "update"
                elif line.startswith("Tasks which are up-to-date"):
                    task_status = "ignore"

                if line.startswith("Task = "):
                    if task_name:
                        yield task_name, task_status, list(split_by_job(block))
                    block = []
                    task_name = re.match("Task = (.*)", line).groups()[0]
                    continue
                if line:
                    block.append(line)
            if task_name:
                yield task_name, task_status, list(split_by_job(block))

        # create connection
        try:
            connection = pika.BlockingConnection(pika.ConnectionParameters(
                host=host))
            self.connected = True
        except pika.exceptions.AMQPConnectionError:
            self.connected = False
            return

        self.channel = connection.channel()
        self.channel.exchange_declare(
            exchange=self.exchange,
            type='topic')

        # populate with initial messages
        for task_name, task_status, jobs in split_by_task(ruffus_text):
            if task_name.startswith("(mkdir"):
                continue

            to_run = 0
            for job_name, job_status, job_message in jobs:
                self.jobs[job_name] = (task_name, job_name)
                if job_status == "update":
                    to_run += 1

            self.tasks[task_name] = [task_status, len(jobs),
                                     len(jobs) - to_run]
            self.send_task(task_name)

    def send_task(self, task_name):
        '''send task status.'''

        if not self.connected:
            return

        task_status, task_total, task_completed = self.tasks[task_name]

        data = {}
        data['created_at'] = time.time()
        data['pipeline'] = self.pipeline_name
        data['task_name'] = task_name
        data['task_status'] = task_status
        data['task_total'] = task_total
        data['task_completed'] = task_completed

        key = "%s.%s.%s" % (self.project_name, self.pipeline_name, task_name)

        try:
            self.channel.basic_publish(exchange=self.exchange,
                                       routing_key=key,
                                       body=json.dumps(data))
        except pika.exceptions.ConnectionClosed:
            E.warn("could not send message - connection closed")
        except Exception as e:
            E.warn("could not send message: %s" % str(e))

    def send_error(self, task_name, job, error=None, msg=None):

        if not self.connected:
            return

        try:
            task_status, task_total, task_completed = self.tasks[task_name]
        except KeyError:
            E.warn("could not get task information for %s, no message sent" %
                   task_name)
            return

        data = {}
        data['created_at'] = time.time()
        data['pipeline'] = self.pipeline_name
        data['task_name'] = task_name
        data['task_status'] = 'failed'
        data['task_total'] = task_total
        data['task_completed'] = task_completed

        key = "%s.%s.%s" % (self.project_name, self.pipeline_name, task_name)

        try:
            self.channel.basic_publish(exchange=self.exchange,
                                       routing_key=key,
                                       body=json.dumps(data))
        except pika.exceptions.ConnectionClosed:
            E.warn("could not send message - connection closed")
        except Exception as e:
            E.warn("could not send message: %s" % str(e))

    def filter(self, record):

        if not self.connected:
            return True

        # filter ruffus logging messages
        if record.filename.endswith("task.py"):
            try:
                before, task_name = record.msg.split(" = ")
            except ValueError:
                return True

            # ignore the mkdir, etc tasks
            if task_name not in self.tasks:
                return True

            if before == "Task enters queue":
                self.tasks[task_name][0] = "running"
            elif before == "Completed Task":
                self.tasks[task_name][0] = "completed"
            elif before == "Uptodate Task":
                self.tasks[task_name][0] = "uptodate"
            else:
                return True

            # send new task status out
            self.send_task(task_name)

        return True


USAGE = '''
usage: %prog [OPTIONS] [CMD] [target]

Execute pipeline %prog.

Commands can be any of the following

make <target>
   run all tasks required to build *target*

show <target>
   show tasks required to build *target* without executing them

plot <target>
   plot image (using inkscape) of pipeline state for *target*

debug <target> [args]
   debug a method using the supplied arguments. The method <target>
   in the pipeline is run without checking any dependencies.

config
   write new configuration files pipeline.ini, sphinxreport.ini and conf.py
   with default values

dump
   write pipeline configuration to stdout

touch
   touch files only, do not run

regenerate
   regenerate the ruffus checkpoint file

check
   check if requirements (external tool dependencies) are satisfied.

clone <source>
   create a clone of a pipeline in <source> in the current
   directory. The cloning process aims to use soft linking to files
   (not directories) as much as possible.  Time stamps are
   preserved. Cloning is useful if a pipeline needs to be re-run from
   a certain point but the original pipeline should be preserved.

'''


def main(args=sys.argv):
    """command line control function for a pipeline.

    This method defines command line options for the pipeline and
    updates the global configuration dictionary correspondingly.

    It then provides a command parser to execute particular tasks
    using the ruffus pipeline control functions. See the generated
    command line help for usage.

    To use it, add::

        import CGAT.Pipeline as P

        if __name__ == "__main__":
            sys.exit(P.main(sys.argv))

    to your pipeline script.

    Arguments
    ---------
    args : list
        List of command line arguments.

    """

    global GLOBAL_OPTIONS
    global GLOBAL_ARGS

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=USAGE)

    parser.add_option("--pipeline-action", dest="pipeline_action",
                      type="choice",
                      choices=(
                          "make", "show", "plot", "dump", "config", "clone",
                          "check", "regenerate"),
                      help="action to take [default=%default].")

    parser.add_option("--pipeline-format", dest="pipeline_format",
                      type="choice",
                      choices=("dot", "jpg", "svg", "ps", "png"),
                      help="pipeline format [default=%default].")

    parser.add_option("-n", "--dry-run", dest="dry_run",
                      action="store_true",
                      help="perform a dry run (do not execute any shell "
                      "commands) [default=%default].")

    parser.add_option("-f", "--force-output", dest="force",
                      action="store_true",
                      help="force running the pipeline even if there "
                      "are uncommited changes "
                      "in the repository [default=%default].")

    parser.add_option("-p", "--multiprocess", dest="multiprocess", type="int",
                      help="number of parallel processes to use on "
                      "submit host "
                      "(different from number of jobs to use for "
                      "cluster jobs) "
                      "[default=%default].")

    parser.add_option("-e", "--exceptions", dest="log_exceptions",
                      action="store_true",
                      help="echo exceptions immediately as they occur "
                      "[default=%default].")

    parser.add_option("-i", "--terminate", dest="terminate",
                      action="store_true",
                      help="terminate immediately at the first exception "
                      "[default=%default].")

    parser.add_option("-d", "--debug", dest="debug",
                      action="store_true",
                      help="output debugging information on console, "
                      "and not the logfile "
                      "[default=%default].")

    parser.add_option("-s", "--set", dest="variables_to_set",
                      type="string", action="append",
                      help="explicitely set paramater values "
                      "[default=%default].")

    parser.add_option("-c", "--checksums", dest="ruffus_checksums_level",
                      type="int",
                      help="set the level of ruffus checksums"
                      "[default=%default].")

    parser.add_option("-t", "--is-test", dest="is_test",
                      action="store_true",
                      help="this is a test run"
                      "[default=%default].")

    parser.add_option("--rabbitmq-exchange", dest="rabbitmq_exchange",
                      type="string",
                      help="RabbitMQ exchange to send log messages to "
                      "[default=%default].")

    parser.add_option("--rabbitmq-host", dest="rabbitmq_host",
                      type="string",
                      help="RabbitMQ host to send log messages to "
                      "[default=%default].")

    parser.set_defaults(
        pipeline_action=None,
        pipeline_format="svg",
        pipeline_targets=[],
        multiprocess=40,
        logfile="pipeline.log",
        dry_run=False,
        force=False,
        log_exceptions=False,
        exceptions_terminate_immediately=False,
        debug=False,
        variables_to_set=[],
        is_test=False,
        ruffus_checksums_level=0,
        rabbitmq_host="saruman",
        rabbitmq_exchange="ruffus_pipelines")

    (options, args) = E.Start(parser,
                              add_cluster_options=True)

    GLOBAL_OPTIONS, GLOBAL_ARGS = options, args
    E.info("Started in: %s" % PARAMS.get("workingdir"))
    # At this point, the PARAMS dictionary has already been
    # built. It now needs to be updated with selected command
    # line options as these should always take precedence over
    # configuration files.

    PARAMS["dryrun"] = options.dry_run
    if options.cluster_queue is not None:
        PARAMS["cluster_queue"] = options.cluster_queue
    if options.cluster_priority is not None:
        PARAMS["cluster_priority"] = options.cluster_priority
    if options.cluster_num_jobs is not None:
        PARAMS["cluster_num_jobs"] = options.cluster_num_jobs
    if options.cluster_options is not None:
        PARAMS["cluster_options"] = options.cluster_options
    if options.cluster_parallel_environment is not None:
        PARAMS["cluster_parallel_environment"] =\
            options.cluster_parallel_environment

    PARAMS["ruffus_checksums_level"] = options.ruffus_checksums_level

    for variables in options.variables_to_set:
        variable, value = variables.split("=")
        PARAMS[variable.strip()] = IOTools.str2val(value.strip())

    version = None

    try:
        # this is for backwards compatibility
        # get mercurial version
        repo = hgapi.Repo(PARAMS["pipeline_scriptsdir"])
        version = repo.hg_id()

        status = repo.hg_status()
        if status["M"] or status["A"]:
            if not options.force:
                raise ValueError(
                    ("uncommitted change in code "
                     "repository at '%s'. Either commit or "
                     "use --force-output") % PARAMS["pipeline_scriptsdir"])
            else:
                E.warn("uncommitted changes in code repository - ignored ")
        version = version[:-1]
    except:
        # try git:
        try:
            stdout, stderr = execute(
                "git rev-parse HEAD", cwd=PARAMS["pipeline_scriptsdir"])
        except:
            stdout = "NA"
        version = stdout

    if args:
        options.pipeline_action = args[0]
        if len(args) > 1:
            options.pipeline_targets.extend(args[1:])

    if options.pipeline_action == "check":
        counter, requirements = Requirements.checkRequirementsFromAllModules()
        for requirement in requirements:
            E.info("\t".join(map(str, requirement)))
        E.info("version check summary: %s" % str(counter))
        E.Stop()
        return

    elif options.pipeline_action == "debug":
        # create the session proxy
        startSession()

        method_name = options.pipeline_targets[0]
        caller = getCaller()
        method = getattr(caller, method_name)
        method(*options.pipeline_targets[1:])

    elif options.pipeline_action in ("make", "show", "svg", "plot",
                                     "touch", "regenerate"):

        # set up extra file logger
        handler = logging.FileHandler(filename=options.logfile,
                                      mode="a")
        handler.setFormatter(
            MultiLineFormatter(
                '%(asctime)s %(levelname)s %(module)s.%(funcName)s.%(lineno)d %(message)s'))
        logger = logging.getLogger()
        logger.addHandler(handler)
        messenger = None

        try:
            if options.pipeline_action == "make":

                # get tasks to be done. This essentially replicates
                # the state information within ruffus.
                stream = StringIO()
                pipeline_printout(
                    stream,
                    options.pipeline_targets,
                    verbose=5,
                    checksum_level=options.ruffus_checksums_level)

                messenger = LoggingFilterRabbitMQ(
                    stream.getvalue(),
                    project_name=getProjectName(),
                    pipeline_name=getPipelineName(),
                    host=options.rabbitmq_host,
                    exchange=options.rabbitmq_exchange)

                logger.addFilter(messenger)

                if not options.without_cluster:
                    global task
                    # use threading instead of multiprocessing in order to
                    # limit the number of concurrent jobs by using the
                    # GIL
                    #
                    # Note that threading might cause problems with rpy.
                    task.Pool = ThreadPool

                    # create the session proxy
                    startSession()

                #
                #   make sure we are not logging at the same time in
                #   different processes
                #
                # session_mutex = manager.Lock()
                E.info(E.GetHeader())
                E.info("code location: %s" % PARAMS["pipeline_scriptsdir"])
                E.info("code version: %s" % version)
                E.info("Working directory is: %s" % PARAMS["workingdir"])

                pipeline_run(
                    options.pipeline_targets,
                    multiprocess=options.multiprocess,
                    logger=logger,
                    verbose=options.loglevel,
                    log_exceptions=options.log_exceptions,
                    exceptions_terminate_immediately=options.exceptions_terminate_immediately,
                    checksum_level=options.ruffus_checksums_level,
                )

                E.info(E.GetFooter())

                closeSession()

            elif options.pipeline_action == "show":
                pipeline_printout(
                    options.stdout,
                    options.pipeline_targets,
                    verbose=options.loglevel,
                    checksum_level=options.ruffus_checksums_level)

            elif options.pipeline_action == "touch":
                pipeline_run(
                    options.pipeline_targets,
                    touch_files_only=True,
                    verbose=options.loglevel,
                    checksum_level=options.ruffus_checksums_level)

            elif options.pipeline_action == "regenerate":
                pipeline_run(
                    options.pipeline_targets,
                    touch_files_only=options.ruffus_checksums_level,
                    verbose=options.loglevel)

            elif options.pipeline_action == "svg":
                pipeline_printout_graph(
                    options.stdout,
                    options.pipeline_format,
                    options.pipeline_targets,
                    checksum_level=options.ruffus_checksums_level)

            elif options.pipeline_action == "plot":
                outf, filename = tempfile.mkstemp()
                pipeline_printout_graph(
                    os.fdopen(outf, "w"),
                    options.pipeline_format,
                    options.pipeline_targets,
                    checksum_level=options.ruffus_checksums_level)
                execute("inkscape %s" % filename)
                os.unlink(filename)

        except ruffus_exceptions.RethrownJobError, value:

            if not options.debug:
                E.error("%i tasks with errors, please see summary below:" %
                        len(value.args))
                for idx, e in enumerate(value.args):
                    task, job, error, msg, traceback = e

                    if task is None:
                        # this seems to be errors originating within ruffus
                        # such as a missing dependency
                        # msg then contains a RethrownJobJerror
                        msg = str(msg)
                        pass
                    else:
                        task = re.sub("__main__.", "", task)
                        job = re.sub("\s", "", job)

                    if messenger:
                        messenger.send_error(task, job, error, msg)

                    # display only single line messages
                    if len([x for x in msg.split("\n") if x != ""]) > 1:
                        msg = ""

                    E.error("%i: Task=%s Error=%s %s: %s" %
                            (idx, task, error, job, msg))

                E.error("full traceback is in %s" % options.logfile)

                # write full traceback to log file only by removing the stdout
                # handler
                lhStdout = logger.handlers[0]
                logger.removeHandler(lhStdout)
                logger.error("start of error messages")
                logger.error(value)
                logger.error("end of error messages")
                logger.addHandler(lhStdout)

                # raise error
                raise ValueError(
                    "pipeline failed with %i errors" % len(value.args))
            else:
                raise

    elif options.pipeline_action == "dump":
        # convert to normal dictionary (not defaultdict) for parsing purposes
        # do not change this format below as it is exec'd in peekParameters()
        print "dump = %s" % str(dict(PARAMS))

    elif options.pipeline_action == "config":
        f = sys._getframe(1)
        caller = inspect.getargvalues(f).locals["__file__"]
        prefix = os.path.splitext(caller)[0]
        writeConfigFiles(prefix)

    elif options.pipeline_action == "clone":
        clonePipeline(options.pipeline_targets[0])

    else:
        raise ValueError("unknown pipeline action %s" %
                         options.pipeline_action)

    E.Stop()


