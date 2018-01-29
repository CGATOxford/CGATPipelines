'''Pipeline.py - Tools for ruffus pipelines
===========================================

The :mod:`Pipeline` module contains various utility functions for
interfacing CGAT ruffus pipelines with an HPC cluster, uploading data
to databases, providing parameterization, and more.

It is a collection of utility functions covering the topics:

* `Pipeline control`_
* `Logging`_
* `Parameterisation`_
* `Running tasks`_
* `Database upload`_
* `Report building`_

See :doc:`pipelines/pipeline_template` for a pipeline illustrating the
use of this module. See :ref:`PipelineSettingUp` on how to set up a
pipeline.

Pipeline control
----------------

:mod:`Pipeline` provides a :func:`main` function that provides command
line control to a pipeline. To use it, add::

    import CGAT.Pipeline as P
    # ...

    if __name__ == "__main__":
        sys.exit(P.main(sys.argv))

to your pipeline script. Typing::

    python my_pipeline.py --help

will provide the following output:

.. program-output:: python ../CGATPipelines/pipeline_template.py --help

Documentation on using pipelines is at :ref:`PipelineRunning`.

Logging
-------

Logging is set up by :func:`main`. Logging messages will be sent to
the file :file:`pipeline.log` in the current directory.  Additionally,
messages are sent to an RabbitMQ_ message exchange to permit
monitoring of pipeline progress.

Running tasks
-------------

:mod:`Pipeline` provides a :func:`Pipeline.run` method to control
running commandline tools. The :func:`Pipeline.run` method takes care
of distributing these tasks to the cluster. It takes into
consideration command line options such as ``--cluster-queue``. The
command line option ``--local`` will run jobs locally for testing
purposes.

For running Python code that is inside a module in a distributed
function, use the :func:`submit` function. The :func:`execute` method
runs a command locally.

Functions such as :func:`shellquote`, :func:`getCallerLocals`,
:func:`getCaller`, :func:`buildStatement`, :func:`expandStatement`,
:func:`joinStatements` support the parameter interpolation mechanism
used in :mod:`Pipeline`.

Parameterisation
----------------

:mod:`Pipeline` provides hooks for reading pipeline configuration
values from :file:`.ini` files and making them available inside ruffus_
tasks. The fundamental usage is a call to :func:`getParamaters` with
a list of configuration files, typically::

    PARAMS = P.getParameters(
        ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
         "../pipeline.ini",
         "pipeline.ini"])

The :mod:`Pipeline` module defines a global variable :data:`PARAMS`
that provides access the configuration values. To get a handle to
this variable outside a pipeline script, call :func:`getParams`::

    my_cmd = "%(scripts_dir)s/bam2bam.py" % P.getParams()

Functions such as :func:`configToDictionary`, :func:`loadParameters`
:func:`matchParameter`, :func:`substituteParameters` support this
functionality.

Functions such as :func:`asList` and :func:`isTrue` are useful to work
with parameters.

The method :func:`peekParameters` allows one to programmatically read the
parameters of another pipeline.

Temporary files
---------------

Tasks containg multiple steps often require temporary memory storage
locations.  The functions :func:`getTempFilename`, :func:`getTempFile`
and :func:`getTempDir` provide these. These functions are aware of the
temporary storage locations either specified in configuration files or
on the command line and distinguish between the ``private`` locations
that are visible only within a particular compute node, and ``shared``
locations that are visible between compute nodes and typically on a
network mounted location.

Requirements
------------

The methods :func:`checkExecutables`, :func:`checkScripts` and
:func:`checkParameter` check for the presence of executables, scripts
or parameters. These methods are useful to perform pre-run checks
inside a pipeline if a particular requirement is met. But see also the
``check`` commandline command.

Database upload
---------------

To assist with uploading data into a database, :mod:`Pipeline` provides
several utility functions for conveniently uploading data. The :func:`load`
method uploads data in a tab-separated file::

    @transform("*.tsv.gz", suffix(".tsv.gz"), ".load")
    def loadData(infile, outfile):
        P.load(infile, outfile)

The methods :func:`mergeAndLoad` and :func:`concatenateAndLoad` upload
multiple files into same database by combining them first. The method
:func:`createView` creates a table or view derived from other tables
in the database. The function :func:`importFromIterator` uploads
data from a python list or other iterable directly.

The functions :func:`tablequote` and :func:`toTable` translate track
names derived from filenames into names that are suitable for tables.

The method :func:`build_load_statement` can be used to create an
upload command that can be added to command line statements to
directly upload data without storing an intermediate file.

The method :func:`connect` returns a database handle for querying the
database.

Report building
---------------

Once built, a report can be published by copying it to the publicly
visible directories on the CGAT systems. At the same time, references
to files on CGAT systems need to be replaced with links through the
public web interface. The function :func:`publish_report` implements
this functionality.

The function :meth:`publish_tracks` builds a UCSC track hub and moves
it into the appropriate CGAT download directories. The method
:func:`publish_notebooks` builds and exports the ipython_ notebooks
related to project.

For these methods to work, the code assumes a certain directory
layout. The method :func:`isCGAT` checks if the code is executed within the
CGAT systems. The functions :func:`getProjectDirectories`,
:func:`getPipelineName`, :func:`getProjectId`, :func:`getProjectName`
provide information about the pipeline executed and the project context.

.. note::

   The methods above are CGAT specific and require a particalur layout.

Package layout
--------------

The module is arranged as a python package with several submodules. Functions
within a submodule to be exported are all imported to the namespace of
:mod:`Pipeline`.

.. toctree::

   Pipeline/Control
   Pipeline/Database
   Pipeline/Execution
   Pipeline/Files
   Pipeline/Local
   Pipeline/Parameters
   Pipeline/Utils

Reference
---------

'''
import os
import sys
import pickle

# import submodules into namespace
from CGATPipelines.Pipeline.Control import *
from CGATPipelines.Pipeline.Database import *
from CGATPipelines.Pipeline.Local import *
from CGATPipelines.Pipeline.Files import *
from CGATPipelines.Pipeline.Cluster import *
from CGATPipelines.Pipeline.Execution import *
from CGATPipelines.Pipeline.Utils import *
from CGATPipelines.Pipeline.Parameters import *


from CGAT import Experiment as E

# import into namespace for backwards compatibility
from CGAT.IOTools import cloneFile as clone
from CGAT.IOTools import touchFile as touch
from CGAT.IOTools import snip as snip

# import submodules
from . import Local as Local
from . import Execution as Execution
from . import Control as Control
from . import Database as Database
from . import Files as Files
from . import Parameters as Parameters

# broadcast parameters and config object, take from
# Parameters.py
PARAMS = Parameters.PARAMS
CONFIG = Parameters.CONFIG

# and drop PARAMS/CONFIG variables into the submodules
Local.CONFIG = CONFIG
Local.PARAMS = PARAMS
Control.PARAMS = PARAMS
Execution.PARAMS = PARAMS
Files.PARAMS = PARAMS

# set working directory at process launch to prevent repeated calls to
# os.getcwd failing if network is busy
PARAMS["workingdir"] = os.getcwd()


def run_report(clean=True,
               with_pipeline_status=True,
               pipeline_status_format="svg"):
    '''run CGATreport.

    This will also run ruffus to create an svg image of the pipeline
    status unless *with_pipeline_status* is set to False. The image
    will be saved into the export directory.

    '''

    if with_pipeline_status:
        targetdir = PARAMS["exportdir"]
        if not os.path.exists(targetdir):
            os.mkdir(targetdir)

        pipeline_printout_graph(
            os.path.join(
                targetdir,
                "pipeline.%s" % pipeline_status_format),
            pipeline_status_format,
            ["full"],
            checksum_level=PARAMS["ruffus_checksums_level"]
        )

    dirname, basename = os.path.split(getCaller().__file__)

    report_engine = PARAMS.get("report_engine", "cgatreport")
    assert report_engine in ('sphinxreport', 'cgatreport')

    docdir = os.path.join(dirname, "pipeline_docs", snip(basename, ".py"))
    themedir = os.path.join(dirname, "pipeline_docs", "themes")
    relpath = os.path.relpath(docdir)
    trackerdir = os.path.join(docdir, "trackers")

    # warning: memory gets multiplied by threads, so set it not too
    # high
    job_memory = PARAMS["report_memory"]
                 #"1G" # This causes problems in outside HPCs

    job_threads = PARAMS["report_threads"]

    # use a fake X display in order to avoid windows popping up
    # from R plots.
    xvfb_command = IOTools.which("xvfb-run")

    # permit multiple servers using -a option
    if xvfb_command:
        xvfb_command += " -a "
    else:
        xvfb_command = ""

    # if there is no DISPLAY variable set, xvfb runs, but
    # exits with error when killing process. Thus, ignore return
    # value.
    # print os.getenv("DISPLAY"), "command=", xvfb_command
    if not os.getenv("DISPLAY"):
        erase_return = "|| true"
    else:
        erase_return = ""

    # in the current version, xvfb always returns with an error, thus
    # ignore these.
    erase_return = "|| true"

    if clean:
        clean = """rm -rf report _cache _static;"""
    else:
        clean = ""

    # with sphinx >1.3.1 the PYTHONPATH needs to be set explicitely as
    # the virtual environment seems to be stripped. It is thus set to
    # the contents of the current sys.path
    syspath = ":".join(sys.path)

    statement = '''
    %(clean)s
    (export SPHINX_DOCSDIR=%(docdir)s;
    export SPHINX_THEMEDIR=%(themedir)s;
    export PYTHONPATH=%(syspath)s;
    %(xvfb_command)s
    %(report_engine)s-build
    --num-jobs=%(report_threads)s
    sphinx-build
    -b html
    -d %(report_doctrees)s
    -c .
    -j %(report_threads)s
    %(docdir)s %(report_html)s
    >& report.log %(erase_return)s )
    '''

    run()

    E.info('the report is available at %s' % os.path.abspath(
        os.path.join(PARAMS['report_html'], "contents.html")))


def publish_notebooks():
    '''publish report into web directory.'''

    dirs = getProjectDirectories()

    notebookdir = dirs['notebookdir']
    exportdir = dirs['exportdir']
    exportnotebookdir = os.path.join(exportdir, "notebooks")

    if not os.path.exists(exportnotebookdir):
        os.makedirs(exportnotebookdir)

    statement = '''
    cd %(exportnotebookdir)s;
    ipython nbconvert
    %(notebookdir)s/*.ipynb
    --to html
    ''' % locals()

    E.run(statement)

__all__ = [
    # backwards incompatibility
    "clone",
    "touch",
    "snip",
    # Execution.py
    "run",
    "execute",
    "shellquote",
    "buildStatement",
    "submit",
    "joinStatements",
    "cluster_runnable",
    "run_pickled",
    # Database.py
    "tablequote",
    "toTable",
    "build_load_statement",
    "load",
    "concatenateAndLoad",
    "mergeAndLoad",
    "connect",
    "createView",
    "getDatabaseName",
    "importFromIterator",
    # Utils.py
    "add_doc",
    "isTest",
    "getCallerLocals",
    "getCaller",
    # Control.py
    "main",
    "peekParameters",
    # Files.py
    "getTempFile",
    "getTempDir",
    "getTempFilename",
    "checkScripts",
    "checkExecutables",
    # Local.py
    "run_report",
    "publish_report",
    "publish_notebooks",
    "publish_tracks",
    "getProjectDirectories",
    "getPipelineName",
    "getProjectId",
    "getProjectName",
    "isCGAT",
    # Parameters.py
    "getParameters",
    "loadParameters",
    "matchParameter",
    "substituteParameters",
    "asList",
    "checkParameter",
    "isTrue",
    "configToDictionary",
]
