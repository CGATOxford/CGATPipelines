====================
Using CGAT pipelines
====================

This section provides a tutorial-like introduction to CGAT pipelines.

Introduction
=============

A pipeline takes input data and performs a series of automated steps
(:term:`task`) on it to produce some output data.

Each pipeline is usually coupled with a :term:`CGATReport` document to
summarize and visualize the results.

It really helps if you are familiar with:

   * the unix command line to run and debug the pipeline
   * python_ in order to understand what happens in the pipeline
   * ruffus_ in order to understand the pipeline code
   * sge_ in order to monitor your jobs
   * git_ in order to up-to-date code

.. _PipelineSettingUp:

Setting up a pipeline
======================

**Step 1**: Install the cgat scripts and CGAP Pipelines:

Check that your computing environment is appropriate and follow installation instructions (see :ref:`CGATSetup`).
The directory in which the CGAT code repository is located is the
:term:`source directory`. It will be abbreviated ``<cgat>`` in the
following commands. 

The source directory will contain the pipeline master script named
:file:`CGATPipelines/pipeline_<name>.py`

The default configuration files will be contained in the folder
:file:`CGATPipelines/Pipeline<Name>/`

All our pipelines are written to be lightweight. Therefore, a module file
assoaiated with the pipeline master script, typically named
:file:`CGATPipelines/Pipeline<Name>.py`, is usually where code required to run the tasks
of the pipleine is located. 

**Step 2**: To run a pipeline you will need to create a :term:`working directory`
and enter it. For example::

   mkdir version1
   cd version1/

This is where the pipeline will be executed and files will be generated in this
directory.

**Step 3**: Our pipelines are written to be ran with minimal hard coded
options. Therefore, to run a pipeline an initial configuration file needs to be
generated. A configuration file with all the default values can be obtained by
running::

      cgatflow <name> config

For example, if you wanted to run the QC pipeline you would run::

      cgatflow readqc config


This will create a new :file:`pipeline.ini` file. **YOU MUST EDIT THIS
FILE**. The default values are likely to use the wrong genome or
point to non-existing locations of indices and databases. The
configuration file should be well documented and the format is
simple. The documenation for the `ConfigParser
<http://docs.python.org/library/configparser.html>`_ python module
contains the full specification.

**Step 4**: Add the input files. The required input is specific for each
pipeline in the documentation stringat the; read the pipeline documentation to find out exactly which
files are needed and where they should be put. Commonly, a pipeline
works from input files linked into the :term:`working directory` and
named following pipeline specific conventions.

**Step 5**: You can check if all the external dependencies to tools and
R packages are satisfied by running::

      cgatflow <name> check

See :ref:`ExternalDependencies` for more information.

.. _PipelineRunning:

Running a pipeline
===================

Pipelines are controlled by a single python script called
:file:`pipeline_<name>.py` that lives in the :term:`source
directory`. Command line usage information is available by running::

   cgatflow <name> --help

The basic syntax for ``pipeline_<name>.py`` is::

   cgatflow <name> [workflow options] [workflow arguments]

For example, to run the readqc pipeline you would run the following::

   cgatflow readqc make full

``workflow options`` can be one of the following:

make <task>

   run all tasks required to build :term:`task`

show <task>

   show tasks required to build :term:`task` without executing them

plot <task>

   plot image of workflow (requires `inkscape <http://inkscape.org/>`_) of
   pipeline state for :term:`task`

touch <task>

   touch files without running :term:`task` or its pre-requisites. This sets the 
   timestamps for files in :term:`task` and its pre-requisites such that they will 
   seem up-to-date to the pipeline.

config

   write a new configuration file :file:`pipeline.ini` with
   default values. An existing configuration file will not be
   overwritten.

clone <srcdir>

   clone a pipeline from :file:`srcdir` into the current
   directory. Cloning attempts to conserve disk space by linking.

In case you are running a long pipeline, make sure you start it
appropriately, for example::

   nice -19 nohup cgatflow <name> make full -v5 -c1

This will keep the pipeline running if you close the terminal.

Additional pipeline options
---------------------------

In addition to running the pipeline with default command line options, running a
pipeline with --help will allow you to see additional options for ``workflow arguments``
when running the pipelines. These will modify the way the pipeline in ran.

`- -local`

    This option allows the pipeline to run locally.

`- -input-validation`

    This option will check the pipeline.ini file for missing values before the
    pipeline starts.

`- -debug`

    Add debugging information to the console and not the logfile

`- -dry-run`

    Perform a dry run of the pipeline (do not execute shell commands)

`- -exceptions`

    Echo exceptions immidietly as they occur.

`-c - -checksums`

    Set the level of ruffus checksums.

Building pipeline reports
================================

Some of the pipelines are associated with an automated report
generator to display summary information as a set of nicely formatted
html pages. 

Currently in CGAT we have 4 different types of report generation.

   * MultiQC report
   * R markdown
   * IPython notebook
   * CGAT reports

To determine which type of reporting is implimented for each pipeline, refer to
the specific pipeline documentation at the beginning of the script.

Reports are generated using the following command once a workflow has completed::

    cgatflow <name> make build_report

MultiQC report
--------------

MultiQC is a python framework for automating reporting and we have imliemnted it in the
majority of our workflows to generate QC stats for frequently used tools (mostly in our
generic workflows). 


R markdown
----------
R markdown report generation is very useful for generating bespoke reports that require user
defined reporting. We have implimented this in our bamstats workflow.

Jupyter notebook
----------------
Jupyter notebook is a second approach that we use to produce bespoke reports. An example is
also implimented in our bamstats workflow.

CGAT Reports
------------
CGAT reports in an in house reporting tool that we have used in the majority of our pipelines.
However, we are depricating its use in the future and most reports will be replaced with either
MultiQC, Rmarkdown or Jupyter reports.

To run CGAT reports:

In order to build the documentation, drop the appropriate
:file:`conf.py` and :file:`cgatreport.ini` configuration files into
the :term:`working directory` and run the pipeline command::

   nice -19 cgatflow <name> make build_report

This will create the report from scratch in the current directory. The
report can be viewed opening the file
:file:`<work>/report/html/contents.html` in your browser.

CGATReport is powerful and can take its time on large projects that
need to generate a multitude of plots and tables. In order to speed up
this process, there are some advanced features that CGATReport offers:

   * caching of results
   * multiprocessing
   * incremental builds
   * separate build directory

Please see the CGATReport_ documentation for more information.


Troubleshooting
===============

Many things can go wrong while running the pipeline. Look out for

   * bad input format. The pipeline does not perform sanity checks on
       the input format.  If the input is bad, you might see wrong or
       missing results or an error message.
   * pipeline disruptions. Problems with the cluster, the file system
       or the controlling terminal might all cause the pipeline to
       abort.
   * bugs. The pipeline makes many implicit assumptions about the
       input files and the programs it runs. If program versions
       change or inputs change, the pipeline might not be able to deal
       with it.  The result will be wrong or missing results or an
       error message.

If the pipeline aborts, locate the step that caused the error by
reading the logfiles and the error messages on stderr
(:file:`nohup.out`). See if you can understand the error and guess the
likely problem (new program versions, badly formatted input, ...). If
you are able to fix the error, remove the output files of the step in
which the error occured and restart the pipeline. Processing should
resume at the appropriate point.

.. note:: 

   Look out for upstream errors. For example, the pipeline might build
   a geneset filtering by a certain set of contigs. If the contig
   names do not match, the geneset will be empty, but the geneset
   building step might conclude successfully. However, you might get
   an error in any of the downstream steps complaining that the gene
   set is empty. To fix this, fix the error and delete the files
   created by the geneset building step and not just the step that
   threw the error.

Common pipeline errors
----------------------

One of the most common errors when runnig the pipeline is::

    GLOBAL_SESSION = drmaa.Session()
    NameError: name 'drmaa' is not defined

This error occurrs because you are not connected to the cluster. Alternatively
you can run the pipleine in local mode by adding --local as a command line option.

Updating to the latest code version
-----------------------------------

To get the latest bugfixes, go into the :term:`source directory` and type::

   git pull

The first command retrieves the latest changes from the master
repository and the second command updates your local version with
these changes.

.. _PipelineReporting:


 
