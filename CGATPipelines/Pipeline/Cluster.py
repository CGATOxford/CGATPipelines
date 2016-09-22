##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2016 Stephen Sansom
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
'''Cluster.py - cluster utility functions for ruffus pipelines
==============================================================

This module abstracts the DRMAA native specification and provides
convenience functions for running Drmaa jobs.


Reference
---------

'''

import re
import os
import stat
import time
import CGAT.Experiment as E

try:
    import drmaa
    HAS_DRMAA = True
except RuntimeError:
    HAS_DRMAA = False


def setupDrmaaJobTemplate(drmaa_session, options, job_name, job_memory):
    '''Sets up a Drmma job template. Currently SGE and SLURM are
       supported'''

    if not job_memory:
        raise ValueError("Job memory must be specified when running"
                         "DRMAA jobs")

    jt = drmaa_session.createJobTemplate()
    jt.workingDirectory = options["workingdir"]
    jt.jobEnvironment = {'BASH_ENV': '~/.bashrc'}
    jt.args = []
    if not re.match("[a-zA-Z]", job_name[0]):
        job_name = "_" + job_name

    # queue manager specific configuration options
    queue_manager = options["cluster_queue_manager"]

    if queue_manager.lower() == "sge":

        # see: ? cannot find documentation on the SGE native spec

        spec = ["-V",
                "-N %s" % job_name]

        if options["cluster_priority"]:
            spec.append("-p %(cluster_priority)i")

        if options["cluster_options"]:
            spec.append("%(cluster_options)s")

        if not options["cluster_memory_resource"]:
            raise ValueError("The cluster memory resource must be specified")

        for resource in options["cluster_memory_resource"].split(","):
            spec.append("-l %s=%s" % (resource, job_memory))

        # if process has multiple threads, use a parallel environment
        multithread = 'job_threads' in options and options['job_threads'] > 1
        if multithread:
            spec.append(
                "-pe %(cluster_parallel_environment)s %(job_threads)i -R y")

        if "cluster_pe_queue" in options and multithread:
            spec.append(
                "-q %(cluster_pe_queue)s")
        elif options['cluster_queue'] != "NONE":
            spec.append("-q %(cluster_queue)s")

    elif queue_manager.lower() == "slurm":

        # SLURM DOCS:
        # http://apps.man.poznan.pl/trac/slurm-drmaa
        # https://computing.llnl.gov/linux/slurm/cons_res_share.html
        #
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
        # -V: SLURM automatically passess the environment variables
        # -p: does not appear to be part of the slurm drmaa native spec
        #
        # TODO: add "--account" (not sure the best way to fill param).

        spec = ["-J %s" % job_name]

        if options["cluster_options"]:
            spec.append("%(cluster_options)s")

        if 'job_threads' in options:
            job_threads = options["job_threads"]
        else:
            job_threads = 1  # probably should come from a config option

        spec.append("--cpus-per-task=%s" % job_threads)

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

        spec.append("--mem-per-cpu=%s" % job_memory_per_cpu)

        # set the partition to use (equivalent of SGE queue)
        spec.append("--partition=%(cluster_queue)s")

    elif queue_manager.lower() == "torque":

        # PBS Torque native specifictation:
        # http://apps.man.poznan.pl/trac/pbs-drmaa

        spec = ["-N %s" % job_name,
                "-l mem=%s" % job_memory, ]

        if options["cluster_options"]:
            spec.append("%(cluster_options)s")

        # There is no equivalent to sge -V option for pbs-drmaa
        # recreating this...
        jt.jobEnvironment = os.environ
        jt.jobEnvironment.update({'BASH_ENV': os.path.join(os.environ['HOME'],
                                                           '.bashrc')})

    else:
        raise ValueError("Queue manager %s not supported" % queue_manager)

    jt.nativeSpecification = " ".join(spec) % options

    # keep stdout and stderr separate
    jt.joinFiles = False

    return jt


def setDrmaaJobPaths(job_template, job_path):
    '''Adds the job_path, stdout_path and stderr_paths
       to the job_template.
    '''
    job_path = os.path.abspath(job_path)

    os.chmod(job_path, stat.S_IRWXG | stat.S_IRWXU)

    stdout_path = job_path + ".stdout"
    stderr_path = job_path + ".stderr"

    job_template.remoteCommand = job_path
    job_template.outputPath = ":" + stdout_path
    job_template.errorPath = ":" + stderr_path

    return job_template, stdout_path, stderr_path


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

    if ignore_pipe_errors:
        return statement
    else:
        return " ".join((_exec_prefix, statement, _exec_suffix))


def collectSingleJobFromCluster(session, job_id,
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
