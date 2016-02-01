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
"""Execution.py - Job control for ruffus pipelines
=========================================================

Session
-------

This module manages a DRMAA session. :func:`startSession`
starts a session and :func:`closeSession` closes it.

Reference
---------

"""

import importlib
import os
import pickle
import pipes
import re
import stat
import subprocess
import sys
import time

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from CGAT.IOTools import snip as snip

# talking to a cluster
try:
    import drmaa
    HAS_DRMAA = True
except RuntimeError:
    HAS_DRMAA = False

# Set from Pipeline.py
PARAMS = {}

from CGATPipelines.Pipeline.Utils import getCallerLocals
from CGATPipelines.Pipeline.Parameters import substituteParameters
from CGATPipelines.Pipeline.Files import getTempFilename, getTempFile

# global drmaa session
GLOBAL_SESSION = None


def _pickle_args(args, kwargs):
    ''' Pickle a set of function arguments. Removes any kwargs that are
    arguements to submit first. Returns a tuple, the first member of which
    is the key word arguements to submit, the second is a file name
    with the picked call arguements '''

    use_args = ["to_cluster",
                "logfile",
                "job_options",
                "job_queue",
                "job_threads",
                "job_memory"]

    submit_args = {}

    for arg in use_args:
        if arg in kwargs:
            submit_args[arg] = kwargs[arg]
            del kwargs[arg]

    args_file = getTempFilename(shared=True)
    pickle.dump([args, kwargs], open(args_file, "wb"))
    return (submit_args, args_file)


def startSession():
    """start and initialize the global DRMAA session."""

    global GLOBAL_SESSION
    GLOBAL_SESSION = drmaa.Session()
    GLOBAL_SESSION.initialize()
    return GLOBAL_SESSION


def closeSession():
    """close the global DRMAA session."""

    if GLOBAL_SESSION is not None:
        GLOBAL_SESSION.exit()


def shellquote(statement):
    '''shell quote a string to be used as a function argument.

    from http://stackoverflow.com/questions/967443/python-module-to-shellquote-unshellquote
    '''
    _quote_pos = re.compile('(?=[^-0-9a-zA-Z_./\n])')

    if statement:
        return _quote_pos.sub('\\\\', statement).replace('\n', "'\n'")
    else:
        return "''"


def execute(statement, **kwargs):
    '''execute a statement locally.

    This method implements the same parameter interpolation
    as the function :func:`run`.

    Arguments
    ---------
    statement : string
        Command line statement to run.

    Returns
    -------
    stdout : string
        Data sent to standard output by command
    stderr : string
        Data sent to standard error by command
    '''

    if not kwargs:
        kwargs = getCallerLocals()

    kwargs = dict(PARAMS.items() + kwargs.items())

    E.debug("running %s" % (statement % kwargs))

    if "cwd" not in kwargs:
        cwd = PARAMS["workingdir"]
    else:
        cwd = kwargs["cwd"]

    # cleaning up of statement
    # remove new lines and superfluous spaces and tabs
    statement = " ".join(re.sub("\t+", " ", statement).split("\n")).strip()
    if statement.endswith(";"):
        statement = statement[:-1]

    process = subprocess.Popen(statement % kwargs,
                               cwd=cwd,
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    # process.stdin.close()
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "Child was terminated by signal %i: \n"
            "The stderr was: \n%s\n%s\n" %
            (-process.returncode, stderr, statement))

    return stdout, stderr

# Definition of helper functions for job scripts
# detect_pipe_error(): propagate error of programs not at the end of a pipe
# checkpoint(): exit a set of chained commands (via ;) if the previous
# command failed.
_exec_prefix = '''detect_pipe_error_helper()
    {
    while [ "$#" != 0 ] ; do
        # there was an error in at least one program of the pipe
        if [ "$1" != 0 ] ; then return 1 ; fi
        shift 1
    done
    return 0
    }
    detect_pipe_error() {
    detect_pipe_error_helper "${PIPESTATUS[@]}"
    return $?
    }
    checkpoint() {
        detect_pipe_error;
        if [ $? != 0 ]; then exit 1; fi;
    }
    '''

_exec_suffix = "; detect_pipe_error"


def buildStatement(**kwargs):
    '''build a command line statement with paramater interpolation.

    The skeleton of the statement should be defined in kwargs.  The
    method then applies string interpolation using a dictionary built
    from the global configuration dictionary PARAMS, but augmented by
    `kwargs`. The latter takes precedence.

    Arguments
    ---------
    kwargs : dict
        Keyword arguments that are used for parameter interpolation.

    Returns
    -------
    statement : string
        The command line statement with interpolated parameters.

    Raises
    ------
    ValueError
        If ``statement`` is not a key in `kwargs`.

    '''

    if "statement" not in kwargs:
        raise ValueError("'statement' not defined")

    local_params = substituteParameters(**kwargs)

    # build the statement
    try:
        statement = kwargs.get("statement") % local_params
    except KeyError, msg:
        raise KeyError(
            "Error when creating command: could not "
            "find %s in dictionaries" % msg)
    except ValueError, msg:
        raise ValueError("Error when creating command: %s, statement = %s" % (
            msg, kwargs.get("statement")))

    # cleaning up of statement
    # remove new lines and superfluous spaces and tabs
    statement = " ".join(re.sub("\t+", " ", statement).split("\n")).strip()
    if statement.endswith(";"):
        statement = statement[:-1]

    return statement


def expandStatement(statement, ignore_pipe_errors=False):
    '''add generic commands before and after statement.

    The prefixes and suffixes added are defined in :data:`exec_prefix`
    and :data:`exec_suffix`. The main purpose of these prefixs is to
    provide error detection code to detect errors at early steps in a
    series of unix commands within a pipe.

    Arguments
    ---------
    statement : string
        Command line statement to expand
    ignore_pipe_errors : bool
        If False, do not modify statement.

    Returns
    -------
    statement : string
        The expanded statement.

    '''

    if ignore_pipe_errors:
        return statement
    else:
        return " ".join((_exec_prefix, statement, _exec_suffix))


def joinStatements(statements, infile):
    '''join a chain of statements into a single statement.

    Each statement contains an @IN@ or a @OUT@ placeholder or both.
    These will be replaced by the names of successive temporary files.

    In the first statement, @IN@ is replaced with `infile`.

    The last statement should move @IN@ to outfile.

    Arguments
    ---------
    statements : list
        A list of command line statements.
    infile : string
        Filename of the first data set.

    Returns
    -------
    statement : string
        A single command line statement.

    '''

    prefix = getTempFilename()
    pattern = "%s_%%i" % prefix

    result = []
    for x, statement in enumerate(statements):
        if x == 0:
            s = re.sub("@IN@", infile, statement)
        else:
            s = re.sub("@IN@", pattern % x, statement)

        s = re.sub("@OUT@", pattern % (x + 1), s).strip()

        if s.endswith(";"):
            s = s[:-1]
        result.append(s)

    assert prefix != ""
    result.append("rm -f %s*" % prefix)

    result = "; checkpoint ; ".join(result)
    return result


def getStdoutStderr(stdout_path, stderr_path, tries=5):
    '''get stdout/stderr allowing for same lag.

    Try at most *tries* times. If unsuccessfull, throw OSError

    Removes the files once they are read.

    Returns tuple of stdout and stderr.
    '''
    x = tries
    while x >= 0:
        if os.path.exists(stdout_path):
            break
        time.sleep(1)
        x -= 1

    x = tries
    while x >= 0:
        if os.path.exists(stderr_path):
            break
        time.sleep(1)
        x -= 1

    try:
        stdout = open(stdout_path, "r").readlines()
    except IOError, msg:
        E.warn("could not open stdout: %s" % msg)
        stdout = []

    try:
        stderr = open(stderr_path, "r").readlines()
    except IOError, msg:
        E.warn("could not open stdout: %s" % msg)
        stderr = []

    try:
        os.unlink(stdout_path)
        os.unlink(stderr_path)
    except OSError, msg:
        pass

    return stdout, stderr


def _collectSingleJobFromCluster(session, job_id,
                                 statement,
                                 stdout_path, stderr_path,
                                 job_path,
                                 ignore_errors=False):
    '''runs a single job on the cluster.'''
    try:
        retval = session.wait(
            job_id, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    except Exception, msg:
        # ignore message 24 in PBS code 24: drmaa: Job
        # finished but resource usage information and/or
        # termination status could not be provided.":
        if not msg.message.startswith("code 24"):
            raise
        retval = None

    stdout, stderr = getStdoutStderr(stdout_path, stderr_path)

    if retval and retval.exitStatus != 0 and not ignore_errors:
        raise OSError(
            "---------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was: \n%s\n%s\n"
            "-----------------------------------------" %
            (retval.exitStatus,
             "".join(stderr), statement))

    try:
        os.unlink(job_path)
    except OSError:
        E.warn(
            ("temporary job file %s not present for "
             "clean-up - ignored") % job_path)


def run(**kwargs):
    """run a command line statement.

    The method runs a single or multiple statements on the cluster
    using drmaa. The cluster is bypassed if:

        * ``to_cluster`` is set to None in the context of the
          calling function.

        * ``--local`` has been specified on the command line
          and the option ``without_cluster`` has been set as
          a result.

        * no libdrmaa is present

        * the global session is not initialized (GLOBAL_SESSION is
          None)

    To decide which statement to run, the method works by examining
    the context of the calling function for a variable called
    ``statement`` or ``statements``.

    If ``statements`` is defined, multiple job scripts are created and
    sent to the cluster. If ``statement`` is defined, a single job
    script is created and sent to the cluster. Additionally, if
    ``job_array`` is defined, the single statement will be submitted
    as an array job.

    Troubleshooting:

       1. DRMAA creates sessions and their is a limited number
          of sessions available. If there are two many or sessions
          become not available after failed jobs, use ``qconf -secl``
          to list sessions and ``qconf -kec #`` to delete sessions.

       2. Memory: 1G of free memory can be requested using the job_memory
          variable: ``job_memory = "1G"``
          If there are error messages like "no available queue", then the
          problem could be that a particular complex attribute has
          not been defined (the code should be ``hc`` for ``host:complex``
          and not ``hl`` for ``host:local``. Note that qrsh/qsub directly
          still works.

    """

    # combine options using correct preference
    options = dict(PARAMS.items())
    options.update(getCallerLocals().items())
    options.update(kwargs.items())

    # insert a few legacy synonyms
    options['cluster_options'] = options.get('job_options',
                                             options['cluster_options'])
    options['cluster_queue'] = options.get('job_queue',
                                           options['cluster_queue'])
    options['without_cluster'] = options.get('without_cluster')

    job_memory = None

    if 'job_memory' in options:
        job_memory = options['job_memory']

    elif "mem_free" in options["cluster_options"] and \
         PARAMS.get("cluster_memory_resource", False):

        E.warn("use of mem_free in job options is deprecated, please"
               " set job_memory local var instead")

        o = options["cluster_options"]
        x = re.search("-l\s*mem_free\s*=\s*(\S+)", o)
        if x is None:
            raise ValueError(
                "expecting mem_free in '%s'" % o)

        job_memory = x.groups()[0]

        # remove memory spec from job options
        options["cluster_options"] = re.sub(
            "-l\s*mem_free\s*=\s*(\S+)", "", o)
    else:
        job_memory = PARAMS.get("cluster_memory_default", "2G")

    def setupJob(session, options, job_memory, job_name):

        jt = session.createJobTemplate()
        jt.workingDirectory = PARAMS["workingdir"]
        jt.jobEnvironment = {'BASH_ENV': '~/.bashrc'}
        jt.args = []
        if not re.match("[a-zA-Z]", job_name[0]):
            job_name = "_" + job_name

        # queue manager specific configuration options

        queue_manager = PARAMS["queue_manager"]

        if queue_manager.lower() == "sge":

            # see: ? cannot find documentation on the SGE native spec

            spec = ["-V",
                    "-p %(cluster_priority)i",
                    "-N %s" % job_name,
                    "%(cluster_options)s"]

            for resource in PARAMS["cluster_memory_resource"].split(","):
                spec.append("-l %s=%s" % (resource, job_memory))

            # if process has multiple threads, use a parallel environment
            if 'job_threads' in options:
                spec.append(
                    "-pe %(cluster_parallel_environment)s %(job_threads)i -R y")
            if "cluster_pe_queue" in options and 'job_threads' in options:
                    spec.append(
                        "-q %(cluster_pe_queue)s")
            else:
                spec.append("-q %(cluster_queue)s")

        elif queue_manager.lower() == "slurm":

            # SLURM DOCS:
            # http://apps.man.poznan.pl/trac/slurm-drmaa
            # https://computing.llnl.gov/linux/slurm/cons_res_share.html
            #

            if 'job_threads' in options:
                job_threads = options["job_threads"]
            else:
                job_threads = 1  # probably should come from a config option

            # Note the that the specified memory must be per CPU
            # for consistency with the implemented SGE approach
            if job_memory.endswith("G"):
                job_memory_per_cpu = int(job_memory[:-1]) * 1000
            elif job_memory.endswith("M"):
                job_memory_per_cpu = int(job_memory[:-1])
            else:
                raise ValueError('job memory unit not recognised for SLURM, '
                                 'must be either "M" (for Mb) or "G" (for Gb),'
                                 ' e.g. 1G or 1000M for 1 Gigabyte of memory')

            # The SLURM Consumable Resource plugin is required
            # The "CR_CPU_Memory" resource must be specified
            #
            # i.e. in slurm.conf:
            # SelectType=select/cons_res
            # SelectTypeParameters=CR_CPU_Memory
            #
            # * Note that --cpus-per-task will actually refer to cores
            #   with the appropriate Node configuration
            #
            # SLURM-DRMAA DOCS - Note that version 1.2 (SVN) is required
            # http://apps.man.poznan.pl/trac/slurm-drmaa
            #
            # Not implemented:
            # -V: cannot find documentation as to what this flag is for
            # -p: does not appear to be part of the slurm drmaa native spec
            #
            # TODO: add "--account" (not sure the best way to fill param).

        spec = ["-J %s" % job_name,
                "--mem-per-cpu=%s" % job_memory_per_cpu,
                "--cpus-per-task=%s" % job_threads,
                "%(cluster_options)s"]

        else:
            raise ValueError("Queue manager %s not supported" % queue_manager)

        jt.nativeSpecification = " ".join(spec) % options
        # keep stdout and stderr separate
        jt.joinFiles = False

        E.debug("Job spec is: %s" % jt.nativeSpecification)

        return jt

    shellfile = os.path.join(PARAMS["workingdir"], "shell.log")

    pid = os.getpid()
    E.debug('task: pid = %i' % pid)

    # connect to global session
    session = GLOBAL_SESSION
    E.debug('task: pid %i: sge session = %s' % (pid, str(session)))

    ignore_pipe_errors = options.get('ignore_pipe_errors', False)
    ignore_errors = options.get('ignore_errors', False)

    # run on cluster if:
    # * to_cluster is not defined or set to True
    # * command line option without_cluster is set to False
    # * an SGE session is present
    run_on_cluster = ("to_cluster" not in options or
                      options.get("to_cluster")) and \
        not options["without_cluster"] and \
        GLOBAL_SESSION is not None

    # SGE compatible job_name
    job_name = re.sub(
        "[:]", "_",
        os.path.basename(options.get("outfile", "ruffus")))

    def buildJobScript(statement, job_memory, job_name):
        '''build job script from statement.

        returns (name_of_script, stdout_path, stderr_path)
        '''

        tmpfile = getTempFile(dir=PARAMS["workingdir"])
        # disabled: -l -O expand_aliases\n" )
        tmpfile.write("#!/bin/bash\n")
        tmpfile.write(
            'echo "%s : START -> %s" >> %s\n' %
            (job_name, tmpfile.name, shellfile))
        # disabled - problems with quoting
        # tmpfile.write( '''echo 'statement=%s' >> %s\n''' %
        # (shellquote(statement), shellfile) )
        tmpfile.write("set | sed 's/^/%s : /' &>> %s\n" %
                      (job_name, shellfile))
        # module list outputs to stderr, so merge stderr and stdout
        tmpfile.write("module list 2>&1 | sed 's/^/%s: /' &>> %s\n" %
                      (job_name, shellfile))
        tmpfile.write("hostname | sed 's/^/%s: /' &>> %s\n" %
                      (job_name, shellfile))
        tmpfile.write("cat /proc/meminfo | sed 's/^/%s: /' &>> %s\n" %
                      (job_name, shellfile))
        tmpfile.write(
            'echo "%s : END -> %s" >> %s\n' %
            (job_name, tmpfile.name, shellfile))

        # restrict virtual memory
        # Note that there are resources in SGE which could do this directly
        # such as v_hmem.
        # Note that limiting resident set sizes (RSS) with ulimit is not
        # possible in newer kernels.
        tmpfile.write("ulimit -v %i\n" % IOTools.human2bytes(job_memory))

        tmpfile.write(
            expandStatement(
                statement,
                ignore_pipe_errors=ignore_pipe_errors) + "\n")
        tmpfile.close()

        job_path = os.path.abspath(tmpfile.name)
        stdout_path = job_path + ".stdout"
        stderr_path = job_path + ".stderr"

        os.chmod(job_path, stat.S_IRWXG | stat.S_IRWXU)

        return (job_path, stdout_path, stderr_path)

    if run_on_cluster:
        # run multiple jobs
        if options.get("statements"):

            statement_list = []
            for statement in options.get("statements"):
                options["statement"] = statement
                statement_list.append(buildStatement(**options))

            if options.get("dryrun", False):
                return

            jt = setupJob(session, options, job_memory, job_name)

            job_ids, filenames = [], []
            for statement in statement_list:
                E.debug("running statement:\n%s" % statement)

                job_path, stdout_path, stderr_path = buildJobScript(statement,
                                                                    job_memory,
                                                                    job_name)

                jt.remoteCommand = job_path
                jt.outputPath = ":" + stdout_path
                jt.errorPath = ":" + stderr_path

                os.chmod(job_path, stat.S_IRWXG | stat.S_IRWXU)

                job_id = session.runJob(jt)
                job_ids.append(job_id)
                filenames.append((job_path, stdout_path, stderr_path))

                E.debug("job has been submitted with job_id %s" % str(job_id))

            E.debug("waiting for %i jobs to finish " % len(job_ids))
            session.synchronize(job_ids, drmaa.Session.TIMEOUT_WAIT_FOREVER,
                                False)

            # collect and clean up
            for job_id, statement, paths in zip(job_ids, statement_list,
                                                filenames):
                job_path, stdout_path, stderr_path = paths
                _collectSingleJobFromCluster(session, job_id,
                                             statement,
                                             stdout_path,
                                             stderr_path,
                                             job_path,
                                             ignore_errors=ignore_errors)

            session.deleteJobTemplate(jt)

        # run single job on cluster - this can be an array job
        else:

            statement = buildStatement(**options)
            E.debug("running statement:\n%s" % statement)

            if options.get("dryrun", False):
                return

            job_path, stdout_path, stderr_path = buildJobScript(statement,
                                                                job_memory,
                                                                job_name)

            jt = setupJob(session, options, job_memory, job_name)

            jt.remoteCommand = job_path
            # later: allow redirection of stdout and stderr to files;
            # can even be across hosts?
            jt.outputPath = ":" + stdout_path
            jt.errorPath = ":" + stderr_path

            if "job_array" in options and options["job_array"] is not None:
                # run an array job
                start, end, increment = options.get("job_array")
                E.debug("starting an array job: %i-%i,%i" %
                        (start, end, increment))
                # sge works with 1-based, closed intervals
                job_ids = session.runBulkJobs(jt, start + 1, end, increment)
                E.debug("%i array jobs have been submitted as job_id %s" %
                        (len(job_ids), job_ids[0]))
                retval = session.synchronize(
                    job_ids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

                stdout, stderr = getStdoutStderr(stdout_path, stderr_path)

            else:
                # run a single job
                job_id = session.runJob(jt)
                E.debug("job has been submitted with job_id %s" % str(job_id))

                _collectSingleJobFromCluster(session, job_id,
                                             statement,
                                             stdout_path,
                                             stderr_path,
                                             job_path,
                                             ignore_errors=ignore_errors)

            session.deleteJobTemplate(jt)
    else:
        # run job locally on cluster
        statement_list = []
        if options.get("statements"):
            for statement in options.get("statements"):
                options["statement"] = statement
                statement_list.append(buildStatement(**options))
        else:
            statement_list.append(buildStatement(**options))

        if options.get("dryrun", False):
            return

        for statement in statement_list:
            E.debug("running statement:\n%s" % statement)

            # process substitution <() and >() does not
            # work through subprocess directly. Thus,
            # the statement needs to be wrapped in
            # /bin/bash -c '...' in order for bash
            # to interpret the substitution correctly.
            if "<(" in statement or ">(" in statement:
                shell = os.environ.get('SHELL', "/bin/bash")
                if "bash" not in shell:
                    raise ValueError(
                        "require bash for advanced shell syntax: <()")
                # Note: pipes.quote is deprecated in Py3, use shlex.quote
                # (not present in Py2.7).
                statement = pipes.quote(statement)
                statement = "%s -c %s" % (shell, statement)

            process = subprocess.Popen(
                expandStatement(
                    statement,
                    ignore_pipe_errors=ignore_pipe_errors),
                cwd=PARAMS["workingdir"],
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

            # process.stdin.close()
            stdout, stderr = process.communicate()

            if process.returncode != 0 and not ignore_errors:
                raise OSError(
                    "---------------------------------------\n"
                    "Child was terminated by signal %i: \n"
                    "The stderr was: \n%s\n%s\n"
                    "-----------------------------------------" %
                    (-process.returncode, stderr, statement))


def submit(module, function, params=None,
           infiles=None, outfiles=None,
           to_cluster=True,
           logfile=None,
           job_options="",
           job_threads=1,
           job_memory=False):
    '''submit a python *function* as a job to the cluster.

    This method runs the script :file:`run_function` using the
    :func:`run` method in this module thus providing the same
    control options as for command line tools.

    Arguments
    ---------
    module : string
        Module name that contains the function. If `module` is
        not part of the PYTHONPATH, an absolute path can be given.
    function : string
        Name of function to execute
    infiles : string or list
        Filenames of input data
    outfiles : string or list
        Filenames of output data
    logfile : filename
        Logfile to provide to the ``--log`` option
    job_options : string
        String for generic job options for the queuing system
    job_threads : int
        Number of slots (threads/cores/CPU) to use for the task
    job_memory : string
        Amount of memory to reserve for the job.

    '''

    if not job_memory:
        job_memory = PARAMS.get("cluster_memory_default", "2G")

    if type(infiles) in (list, tuple):
        infiles = " ".join(["--input=%s" % x for x in infiles])
    else:
        infiles = "--input=%s" % infiles

    if type(outfiles) in (list, tuple):
        outfiles = " ".join(["--output-section=%s" % x for x in outfiles])
    else:
        outfiles = "--output-section=%s" % outfiles

    if logfile:
        logfile = "--log=%s" % logfile
    else:
        logfile = ""

    if params:
        params = "--params=%s" % ",".join(params)
    else:
        params = ""

    statement = '''python %(pipeline_scriptsdir)s/run_function.py
                          --module=%(module)s
                          --function=%(function)s
                          %(logfile)s
                          %(infiles)s
                          %(outfiles)s
                          %(params)s
                '''
    run()


def cluster_runnable(func):
    '''A dectorator that allows a function to be run on the cluster.

    The decorated function now takes extra arguments. The most important
    is *submit*. If set to true, it will submit the function to the cluster
    via the Pipeline.submit framework. Arguments to the function are
    pickled, so this will only work if arguments are picklable. Other
    arguments to submit are also accepted.

    Note that this allows the unusal combination of *submit* false,
    and *to_cluster* true. This will submit the function as an external
    job, but run it on the local machine.

    Note: all arguments in the decorated function must be passed as
    key-word arguments.
    '''

    # MM: when decorating functions with cluster_runnable, provide
    # them as kwargs, else will throw attribute error

    function_name = func.__name__

    def submit_function(*args, **kwargs):

        if "submit" in kwargs and kwargs["submit"]:
            del kwargs["submit"]
            submit_args, args_file = _pickle_args(args, kwargs)
            module_file = os.path.abspath(
                sys.modules[func.__module__].__file__)
            submit(snip(__file__),
                   "run_pickled",
                   params=[snip(module_file), function_name, args_file],
                   **submit_args)
        else:
            # remove job contral options before running function
            for x in ("submit", "job_options", "job_queue",
                      "job_memory", "job_threads"):
                if x in kwargs:
                    del kwargs[x]
            return func(*args, **kwargs)

    return submit_function


def run_pickled(params):
    ''' run a function whose arguments have been pickled.

    expects that params is [module_name, function_name, arguments_file] '''

    module_name, func_name, args_file = params
    location = os.path.dirname(module_name)
    if location != "":
        sys.path.append(location)

    module_base_name = os.path.basename(module_name)
    E.info("importing module '%s' " % module_base_name)
    E.debug("sys.path is: %s" % sys.path)

    module = importlib.import_module(module_base_name)
    try:
        function = getattr(module, func_name)
    except AttributeError as msg:
        raise AttributeError(msg.message +
                             "unknown function, available functions are: %s" %
                             ",".join([x for x in dir(module)
                                       if not x.startswith("_")]))

    args, kwargs = pickle.load(open(args_file, "rb"))
    E.info("arguments = %s" % str(args))
    E.info("keyword arguments = %s" % str(kwargs))

    function(*args, **kwargs)

    os.unlink(args_file)
