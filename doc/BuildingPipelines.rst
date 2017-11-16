=======================
Building CGAT pipelines
=======================

The best way to build a pipeline is to start from an example. There are several 
pipelines available, see :ref:`cgatpipelines`. To start a new project, use 
:file:`pipeline_quickstart.py`::

   cgatflow quickstart --set-name=test

This will create a report directory and an src directory.

If you navigate to the src directory you will observe that there are two folders
``pipeline_docs/``, ``pipline_test/`` and a ``pipeline_test.py`` pipeline task file.

In order to help with debugging and reading our code, our pipelines are written so that
a pipeline task file contains Ruffus tasks and calls functions in an associated module file,
which contains all of the code to transform and analyse the data.

The module file is not generated during running of the pipeline_testing.py script. Therefore,
if you wish to create a module file, we usually save this file in the following convention,
``PipelineTest.py`` and it can be imported into the main pipeline task file (``pipeline_test.py``).

This section describes how CGAT pipelines can be constructed using the
:mod:`Pipeline` module. The Pipeline.py module contains a variety of
useful functions for pipeline construction.

.. _PipelineOrganization:

Overview
========

Pipelines generally have a similar structure. Pipelines are
implemented as a pipeline script in the :term:`source directory`
called :file:`pipeline_<somename>.py` and a file
:file:`pipeline_<somename>.ini` with default configuration values.

Pipeline input
--------------

Pipelines are executed within a dedicated :term:`working
directory`. They usually require the following files within this
directory:

   * a pipeline configuration file :file:`pipeline.ini`
   * input data files, usually linked in from a data repository

Other files that might be used in a pipeline are:

   * external data files such as genomes that a referred to by they their full path name.
   * conf.py for automated reports.

The pipelines will work from the input files in the :term:`working
directory`, usually identified by their suffix. For example, a
ChIP-Seq pipeline might look for any ``*.fastq.gz`` files in the
directory, run QC on these, map the reads to a genome sequence, call
peaks, do motif analyses, etc.

Pipeline output 
----------------

The pipeline will create files and database tables in the
:term:`working directory`.  When building a pipeline, you can choose
any file/directory layout that suits your needs. Some prefer flat
hierarchies with many files, while others prefer deep directories.

Two directories have a special function and can be used for exporting pipeline results (see PipelinePublishing_):

The :file:`export` directory contains all files that will be referred
to directly in the report or that later should be published by the
pipeline. For example, pdf documents created by the peak caller or
logo images created by a motif tool should go there.

The directory :file:`report` will contain the automatically generated report.

Guidelines
==========

To preserve disk space, please always work use compressed files as
much as possible.  Most data files compress very well, for example
fastq files often compress by a factor of 80% or more: a 10Gb file
will use just 2Gb.

Working with compressed files is straight-forward using unix pipes and
the commands ``gzip``, ``gunzip`` or ``zcat``.

If you require random access to a file, load the file into the
database and index it appropriately. Genomic interval files can be
indexed with tabix to allow random access.

.. _PipelineCommands:

Running commands within tasks
=============================

To run a command line program within a pipeline task, build a
statement and call the :meth:`Pipeline.run` method::

   @files( '*.unsorted', suffix('.unsorted'), '.sorted')
   def sortFile( infile, outfile ):

       statement = '''sort %(infile)s > %(outfile)s'''
       P.run()

On calling the :meth:`Pipeline.run` method, the environment of the
caller is examined for a variable called ``statement``. The variable
is subjected to string substitution from other variables in the local
namespace. In the example above, ``%(infile)s`` and ``%(outfile)s``
are substituted with the values of the variables ``infile`` and
``outfile``, respectively.

The same mechanism also permits setting configuration parameters, for example::

   @files( '*.unsorted', suffix('.unsorted'), '.sorted')
   def sortFile( infile, outfile ):

       statement = '''sort -t %(tmpdir)s %(infile)s > %(outfile)s'''
       P.run()

will automatically substitute the configuration parameter ``tmpdir``
into the command. See ConfigurationValues_ for more on using configuration
parameters.

The pipeline will stop and return an error if the command exits with an error code.

If you chain multiple commands, only the return value of the last
command is used to check for an error. Thus, if an upstream command
fails, it will go unnoticed.  To detect these errors, insert
``&&`` between commands. For example::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted)
   def sortFile( infile, outfile ):

       statement = '''gunzip %(infile)s %(infile)s.tmp &&
		      sort -t %(tmpdir)s %(infile)s.tmp > %(outfile)s &&
		      rm -f %(infile)s.tmp
       P.run()

Of course, the statement aboved could be executed more efficiently
using pipes::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The pipeline inserts code automatically to check for error return
codes if multiple commands are combined in a pipe.

Running commands on the cluster
-------------------------------

In order to run commands on cluster, use ``to_cluster=True``.

To run the command from the previous section on the cluster::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       to_cluster = True
       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The pipeline will automatically create the job submission files,
submit the job to the cluster and wait for its return.

Pipelines will use the command line options ``--cluster-queue``,
``--cluster-priority``, etc. for global job control. For example, to
change the priority when starting the pipeline, use::

   python <pipeline_script.py> --cluster-priority=-20

To set job options specific to a task, you can define additional
variables::

   @files( '*.unsorted.gz', suffix('.unsorted.gz'), '.sorted.gz')
   def sortFile( infile, outfile ):

       to_cluster = True
       job_queue = 'longjobs.q'
       job_priority = -10
       job_options= "-pe dedicated 4 -R y" 
 
       statement = '''gunzip < %(infile)s 
		      | sort -t %(tmpdir)s 
		      | gzip > %(outfile)s'''
       P.run()

The above statement will be run in the queue ``longjobs.q`` at a
priority of ``-10``.  Additionally, it will be executed in the
parallel environment ``dedicated`` with at least 4 cores.

Array jobs can be controlled through the ``job_array`` variable::

   @files( '*.in', suffix('.in'), '.out')
   def myGridTask( infile, outfile ):

       job_array=(0, nsnps, stepsize)
   
       statement = '''grid_task.bash %(infile)s %(outfile)s
          > %(outfile)s.$SGE_TASK_ID 2> %(outfile)s.err.$SGE_TASK_ID
       '''
       P.run()


Note that the :file:`grid_task.bash` file must be grid engine
aware. This means it makes use of the :envvar:`SGE_TASK_ID`,
:envvar:`SGE_TASK_FIRST`, :envvar:`SGE_TASK_LAST` and
:envvar:`SGE_TASK_STEPSIZE` environment variables to select the chunk
of data it wants to work on.

The job submission files are files called `tmp*` in the :term:`working
directory`.  These files will be deleted automatically. However, the
files will remain after aborted runs to be cleaned up manually.

.. _PipelineTracks:

Tracks
======

A pipeline typically processes the data streams from several
experimental data sources. These data streams are usually processed
separately (processing, quality control) and as aggregates. The module
:mod:`PipelineTracks` helps implementing this.

.. _PipelineDatabases:

Databases
=========

Loading data into the database
------------------------------

:mod:`Pipeline.py` offers various tools for working with databases. By
default, it is configured to use an sqlite3 database in the
:term:`working directory` called :file:`csvdb`.

Tab-separated output files can be loaded into a table using the
:meth:`Pipeline.load` function. For example::

   @transform('data_*.tsv.gz', suffix('.tsv.gz'), '.load')
   def loadTables(infile, outfile):
      P.load(infile, outfile)

The task above will load all tables ending with ``tsv.gz`` into the
database Table names are given by the filenames, i.e, the data in
:file:`data_1.tsv.gz` will be loaded into the table :file:`data_1`.

The load mechanism uses the script :file:`csv2db.py` and can be
configured using the configuration options in the ``database`` section
of :file:`pipeline.ini`. Additional options can be given via the
optional *options* argument::

   @transform('data_*.tsv.gz', suffix('.tsv.gz'), '.load')
   def loadTables( infile, outfile ):
      P.load(infile, outfile, "--add-index=gene_id")

In order for the load mechanism to be transparent, it is best avoided
to call the :file:`csv2db.py` script directly. Instead, use the
:meth:`Pipeline.load` function. If the :file:`csv2db.py` needs to
called at the end of a succession of statements, use the
:meth:`Pipeline.build_load_statement` method, for example::

   def loadTranscript2Gene(infile, outfile):
       '''build and load a map of transcript to gene from gtf file
       '''
       load_statement = P.build_load_statement(
	   P.toTable(outfile),
	   options="--add-index=gene_id "
	   "--add-index=transcript_id ")

       statement = '''
       gunzip < %(infile)s
       | python %(scriptsdir)s/gtf2tsv.py --output-map=transcript2gene -v 0
       | %(load_statement)s
       > %(outfile)s'''
       P.run()

See also the variants :meth:`Pipeline.mergeAndLoad` and
`:meth:`Pipeline.concatenateAndLoad` to combine multiple tables and
upload to the database in one go.

Connecting to a database
------------------------

To use data in the database in your tasks, you need to first connect
to the database. The best way to do this is via the connect() method
in Pipeline.py.

The following example illustrates how to use the connection::

    @transform( ... )
    def buildCodingTranscriptSet( infile, outfile ):

	dbh = connect()

	statement = '''SELECT DISTINCT transcript_id FROM transcript_info WHERE transcript_biotype = 'protein_coding' '''
	cc = dbh.cursor()
	transcript_ids = set( [x[0] for x in cc.execute(statement)] )
	...

.. _PipelineReports:

Reports
=======

The :meth:`Pipeline.run_report` method builds or updates reports using
CGATreport_. Usually, a pipeline will simply contain the following::

    @follows( mkdir( "report" ) )
    def build_report():
	'''build report from scratch.'''

	E.info( "starting report build process from scratch" )
	P.run_report( clean = True )

    @follows( mkdir( "report" ) )
    def update_report():
	'''update report.'''

	E.info( "updating report" )
	P.run_report( clean = False )

This will add the two tasks ``build_report`` and ``update_report`` to
the pipeline. The former completely rebuilds a report, while the
latter only updates changed pages. The report will be in the directory
:file:`report`.

Note that report building requires the file :file:`conf.py` in the
:term:`working directory`. This file is read by sphinx_ and can be
used to report building options. By default, the file is a stub
reading in common options from the CGAT code base.

The section :ref:`WritingReports` contains more information.

.. _ConfigurationValues:

Configuration values
====================

Setting up configuration values
--------------------------------

There are different ways to pass on configuration values to pipelines.
Here we explain the priority for all the possible options so you can
choose the best one for your requirements.

The pipeline goes *in order* through different configuration options
to load configuration values and stores them in the :py:data:`PARAMS`
dictionary. This order determines a priority so values read in the first
place can be overwritten by values read in subsequent steps; i.e. values
read lastly have higher priority.

Here is the order in which the configuration values are read:

1. Hard-coded values in :file:`CGATPipelines/Pipeline/Parameters.py`.
2. Parameters stored in :file:`pipeline.ini` files in different locations.
3. Variables declared in the ruffus tasks calling ``P.run()``;
   e.g. ``job_memory=32G``
4. ``cluster_*`` options specified in the command line;
   e.g. ``python pipeline_example.py --cluster-parallel=dedicated make full``

This means that configuration values for the cluster provided as
command-line options will have the highest priority. Therefore::

   python pipeline_example.py --cluster-parallel=dedicated make full

will overwrite any ``cluster_parallel`` configuration values given
in :file:`pipeline.ini` files. Type::

   python pipeline_example.py --help

to check the full list of available command-line options.

You are encouraged to include the following snippet at the beginning
of your pipeline script to setup proper configuration values for
your analyses::

   # load options from the config file
   import CGATPipelines.Pipeline as P
   PARAMS =
     P.getParameters(["%s/pipeline.ini" % os.path.splitext(__file__)[0]])

The method :meth:`Pipeline.getParameters` reads parameters from
the :file:`pipeline.ini` located in the current :term:`working directory`
and updates :py:data:`PARAMS`, a global dictionary of parameter values.
It automatically guesses the type of parameters in the order of ``int()``,
``float()`` or ``str()``. If a configuration variable is empty (``var=``),
it will be set to ``None``.

However, as explained above, there are other :file:`pipeline.ini`
files that are read by the pipeline at start up. In order to get the
priority of them all, you can run::

   python pipeline_example.py printconfig

to see a complete list of :file:`pipeline.ini` files and their priorities.

Configuration values from another pipeline can be added in a separate
namespace::

   PARAMS_ANNOTATIONS = P.peekParameters(
       PARAMS["annotations_dir"],
       "pipeline_annotations.py")

The statement above will load the parameters from a
:mod:`pipeline_annotations` pipeline with :term:`working directory`
``annotations_dir``.

Using configuration values
--------------------------

Configuration values are accessible via the :py:data:`PARAMS`
variable. The :py:data:`PARAMS` variable is a dictionary mapping
configuration parameters to values. Keys are in the format
``section_parameter``. For example, the key ``bowtie_threads`` will
provide the configuration value of::

   [bowtie]
   threads=4

In a script, the value can be accessed via
``PARAMS["bowtie_threads"]``.

Undefined configuration values will throw a :class:`ValueError`. To
test if a configuration variable exists, use::

   if 'bowtie_threads' in PARAMS: pass
      
To test, if it is unset, use::

   if 'bowie_threads' in PARAMS and not PARAMS['botwie_threads']:
      pass

Task specific parameters
------------------------

Task specific parameters can be set by creating a task specific section in
the :file:`pipeline.ini`. The task is identified by the output filename.
For example, given the following task::

   @files( '*.fastq', suffix('.fastq'), '.bam')
   def mapWithBowtie( infile, outfile ):
      ...

and the files :file:`data1.fastq` and :file:`data2.fastq` in the
:term:`working directory`, two output files :file:`data.bam` and
:file:`data2.bam` will be created on executing ``mapWithBowtie``. Both
will use the same parameters. To set parameters specific to the
execution of :file:`data1.fastq`, add the following to
:file:`pipeline.ini`::

   [data1.fastq]
   bowtie_threads=16

This will set the configuration value ``bowtie_threads`` to 16 when
using the command line substitution method in :meth:`Pipeline.run`. To
get an task-specific parameter values in a python task, use::

   @files( '*.fastq', suffix('.fastq'), '.bam')
   def mytask( infile, outfile ):
       MY_PARAMS = P.substituteParameters( locals() )
       
Thus, task specific are implemented generically using the
:meth:`Pipeline.run` mechanism, but pipeline authors need to
explicitely code for track specific parameters.

.. _PipelineDocumentation:

Documentation
=============

Up-to-date and accurate documentation is crucial for writing portable
and maintainable pipelines. To document your pipelines write
documentation as you would for a module.  See
:file:`pipeline_template.py` and other pipelines for an example.

To rebuild all documentation, enter the :file:`doc` directory in the
:term:`source directory` and type::

   cd doc
   python collect.py

This will collect all new scripts to the documentation.

Next, edit the file :file:`contents.rst` and add your pipeline to the
table of pipelines. Finally, type::

   make html

to rebuild the documentation.

Using other pipelines
=====================

You can use the output of other pipelines within your own
pipelines. :mod:`pipeline_annotations` is an example - it provides
often used annotation data sets for an analysis. How to load another
pipelines parameters, connect to its database and write a modular
report have been discussed above.

If you write a pipeline that is likely to be used by others, it is
best to provide an interface.  For example, the
:mod:`pipeline_annotations` pipeline has an interface section that
list all the files that are produced by the pipeline. Other pipelines
can refer to the interface section without having to be aware of the
actual file names::

    filename_cds = os.path.join( PARAMS["annotations_dir"],
             	            PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"] )

Running other pipelines within your pipeline *should* be possible as
well - provided they are within their own separate :term:`working
directory`.

.. _PipelinePublishing:


Checking requisites
===================

Add a ``Requisites`` section to the docstring of the pipeline and any
auxilliary script or module that is called by the pipeline. An example
is below::

    Requirements:

    * tophat >= 2.0.13
    * bowtie2 >= 2.2.3
    * bwa >= 0.7.8
    * gsnap >= 2014-01-21
    * star >= 2.3.0e

These sections allow us to keep track of dependencies for a specific
pipeline and the software collection as a whole.

See :doc:`modules/Requirements` for more information about the
mechanism.

.. _ruffus: http://www.ruffus.org.uk/
.. _sqlite: http://www.sqlite.org/

