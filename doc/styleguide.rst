.. _StyleGuide:

Documentation
=============

Writing doc-strings
-------------------

All pipeline tasks(functions) should be documented through their
doc-string using restructured text.  Pipeline script doc-strings are 
richer and give a more verbose description of the usages, tasks and
targets of the pipeline.  This should also document the different
sections or groups of tasks based on common functionality.

Doc-string lines should be limited to 79 characters (`PEP8 <https://www.python.org/dev/peps/pep-0008/#maximum-line-length>`_).  Sections
are broken by two blank lines

Doc-strings for pipeline scripts
--------------------------------

Each pipeline script doc-string should describe the function(s), usage,
inputs and outputs of the pipeline.  These are broken down into 6 sections::

  * Overview
  * Usage
  * Pipeline output
  * Example
  * Glossary
  * Code (this is generated from the task doc-strings)

Overview
++++++++
This should be a brief overview of the pipeline and its implementation, including
main targets, and optional or ancillary targets under the sub-headers `Principal targets`
and `Optional targets`.

Usage
+++++
This should describe, in detail, the use, configuration and required, or accepted, inputs
of the pipeline.  Sub-headers may include::

  * Configuration
  * Input
  * Option inputs
  * Requirements
  * Additional notes

The `Configuration` sub-section should briefly describe the config file requirements.
The `Input` sub-section should describe any filename format requirements.  For example::

  <Tissue>-<Replicate>-<Condition>.<suffix>

`Requirements` should give a table of third-part tools required for the pipeline to work, 
including the minimum version requirement and the purpose of the tool in the pipeline.
Ideally third-party tool usage should be minimised where possibly to maximise portability
of pipelines, but there's no point in trying to re-invent the wheel either.  For example::

  +---------+------------+------------------------------------------------+
  |*Program*|*Version*   |*Purpose*                                       |
  +---------+------------+------------------------------------------------+
  |bowtie_  |>=0.12.7    |read mapping                                    |
  +---------+------------+------------------------------------------------+
  |tophat_  |>=1.4.0     |read mapping                                    |
  +---------+------------+------------------------------------------------+
  |gsnap_   |>=2012.07.20|read mapping                                    |
  +---------+------------+------------------------------------------------+

Pipeline output
+++++++++++++++
`Pipeline output` describes the major output files from the `Principal targets`
functions.

Example
+++++++
A clear and concise example of the pipeline use should be given in the `Example`
sub-section, including notes of other dependencies not described in the `Requirements`
section.  For example::

  Example
  =======

  Example data is available at
  http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz.  To run
  the example, simply unpack and untar::

    wget http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz
    tar -xvzf pipeline_mapping.tgz
    cd pipeline_mapping
    python <srcdir>/pipeline_mapping.py make full

    .. note::
       For the pipeline to run, install the :doc:`pipeline_annotations` as well.


Glossary
++++++++
`Glossary` gives a set of recognised terms with a short description of their meaning.
For example::

  Glossary
  ========

  .. glossary::

     tophat
       tophat_ - a read mapper to detect splice-junctions

     hisat
       hisat_ - a read mapper for RNASEQ data (basis for tophat3)

     bowtie
        bowtie_ - a read mapper

Code
++++
The `Code` subsection is automatically built from the pipeline task doc-strings.  See below
for the style guide for pipeline tasks.

Writing documentation for pipeline tasks
----------------------------------------

The pipeline tasks(functions) use the NumPy style guide
(`NumPy Style Guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_)
with input Parameters and Returns values.  This sets the minimum standard
for all pipelines.

The doc-string should start with a short description of the task
purpose (1 or 2 lines).  A more verbose description can then follow if
necessary to explain and justify the tasks function and
implementation.  These should be separated by a blank line to denote a
section break.

Task doc-strings always contain the sub-header ``Parameters`` with 
descriptions of the input file, output file and parameters. Each entry in
the ``Parameters`` section has a type, which is the Python data
structure/type to expect, e.g. `list`, `int` or `str`. In addition,
for files, the format should be specified in the description where
appropriate, for example :term:`gtf` or :term:`bed`. If there are
multiple input or output files, each should be described on a separate
line.

Additional ``Parameters`` that are derived from the pipeline config file,
and thus are accessed from the :term:`PARAMS` dictionary should be denoted
in their description by a preceding :term:`PARAMS`.  The description should be
sufficient to give the user the required information to understand the purpose
and implementation of the task.

Optional sub-sections include ``Returns`` where a function returns a
python object. A ``Notes`` sub-section may also be included with
additional information but notes are usually better rendered using the
Sphinx directive ``.. note::`` (see example)


For example::

  def mapMyReads(infiles, outfile):

      '''Map my reads with a short read aligner


      Use MyAligner to map my special short reads against my special reference
      genome.  MyAligner uses MySpecialAlgorithm to align short reads super
      quick against my horrible genome assembly


      Parameters
      ----------
      infiles: str
        the input files for this task - :term:`FASTQ` files generated
	by a short-read sequencer

      mapper_threads: int
        ``PARAMS`` - the number of threads to use for multi-threading mapper

      mapper_memory: str
        ``PARAMS`` - the amount of memory to assign per thread for mapping reads

      outfile: str
        A :term:`BAM` file of reads aligned to a reference genome

	
      .. note::
          Colour space mapping is not implemented

      '''
