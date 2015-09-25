.. _modules:

=======
Modules
=======

Modules fall into two categories, modules that contain functionality
specific to a pipeline, and those that are more widely used.

.. note::

   Functionality in modules should pertain to workflows and
   pipelines. Algorithms, file format parsing and other generic
   functionality should be implemented in scripts or modules in the
   CGAT code collection. See :ref:`Code organization` for more
   information.

Generic modules
===============

Modules in this section contain generic functionality for pipelines
such as for report building, sample handling, etc., and tasks that
occur in multiple pipelines. The most prominent package in this
section is :mod:`Pipeline` that provides command line control,
parameterization and more to pipeline scripts and is organized
into several sub-packages.

.. toctree::
   :maxdepth: 2

   pipelinemodules/Pipeline.rst
   pipelinemodules/PipelineTracks.rst
   pipelinemodules/PipelineUCSC.rst
   pipelinemodules/PipelineGO.rst
   pipelinemodules/PipelineGeneset.rst

Pipeline specific modules
=========================

Modules in this section are pipeline specific. They should be named
after the pipeline they support. Methods in these modules typically
contain the implementation of tasks, tool wrappers and
parsers.

.. toctree::
   :maxdepth: 1

   pipelinemodules/PipelineExome.rst
   pipelinemodules/PipelineIDR.rst
   pipelinemodules/PipelineKEGG.rst
   pipelinemodules/PipelineLncRNA.rst
   pipelinemodules/PipelineMapping.rst
   pipelinemodules/PipelineMappingQC.rst
   pipelinemodules/PipelinePeakcalling.rst
   pipelinemodules/PipelinePreprocess.rst
   pipelinemodules/PipelineReadqc.rst
   pipelinemodules/PipelineRnaseq.rst
   pipelinemodules/PipelineRrbs.rst
   pipelinemodules/PipelineTimeseries.rst
   pipelinemodules/PipelineWindows.rst

Unused modules
==============

The modules below contain tasks that are not in active use, are not
fully developed or are used in very specialized pipelines. They are
a good place to look for existing functionality.

.. toctree::
   :maxdepth: 2

   pipelinemodules/PipelineChipseq.rst
   pipelinemodules/PipelineEnrichment.rst
   pipelinemodules/PipelineIntervalAnnotation.rst
   pipelinemodules/PipelineMedip.rst
   pipelinemodules/PipelineMotifs.rst
   pipelinemodules/PipelineTransfacMatch.rst
   pipelinemodules/PipelineMetagenomeBenchmark.rst
   pipelinemodules/PipelineMetagenomeCommunities.rst
