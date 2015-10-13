=========
Reference
=========

This section contains a reference to the complete contents of the
pipeline collection. It describes the layout of the code repository
and contains a description of the scripts and modules.
collection.

.. _CodeOrganization:

Code organization
=================

The code in the `CGAT Pipeline Collection`_ and the `CGAT Code
Collection`_ is organized on four levels. :term:`pipeline scripts`
and :term:`pipeline modules` are located in the :term:`CGATPipelines`
repository, while :term:`standalone scripts` and :term:`standalone
modules` or located in the `CGAT Code Collection`_.

.. glossary::

   pipeline scripts
   pipeline script

      The :term:`pipeline script`, for example, ``pipeline_mapping.py``
      The pipeline script contains all the functionality that pertains to
      workflow logic. It lists all the tasks that the pipeline
      implements, which tasks should be executed.  It reads the
      configuration parameters and parameterize the workflow and
      individual tasks.  Simple tasks, for example those that call an
      external tool, can be implemented in the :term:`pipeline
      script`. Anything that is more complex should go into
      the :term:`pipeline module`.

   pipeline modules
   pipeline module

      The :term:`pipeline module`, for example,
      ``PipelineMapping.py``. The module associated with a pipeline
      contains implementations of particular tasks. This module should
      not require access to the configuration directory, instead
      parameters should be passed from the pipeline script. The pipeline
      modules should not contain implementations of algorithms, these
      should go into a :term:`standalone script` or :term:`standalone
      module` in the CGAT code collection.

   standalane scripts
   standalone script

      A :term:`standalone script` provides a frontend towards running
      algorithmic procedures on input data. For example, consider the
      task of computing a particular summary statistic from a
      :term:`BAM` formatted file, such as the average mapping quality.
      The looping over the file could be easily achieved with a few
      lines of code::

	  with pysam.AlignmentFile(infile, "rb") as inf:  
	      qualities = [read.mapping_quality for read in inf]
	  mean_quality = sum(qualities) / len(qualities)

      There is a temptation to simply add this piece of code inside a
      :term:`pipeline module` or a :term:`pipeline script`. The best
      place, however, is to add the functionality to a :term:`standalone
      script` as it permits access to the functionality from the command
      line and is thus accessible in many different contexts. Also, a
      :term:`standalone script` can be executed on the cluster. Even
      better, it can be grouped with other functionalities to build
      powerful tools. For example, the above functionality could go into
      a script ``bam2stats.py`` that collects multiple summary statistics
      of :term:`BAM` formatted files.

      Providing a command line front end to even small functionality
      has set-up costs, but the long-term benefits easily offset these.

      The :term:`standalone script` should focus on setting up the
      components and the parameterization of a particular algorithm or
      sequence of computations.  The actual implementation of
      functionality should reside as much as possible in a
      :term:`standalone module`.

   standalone modules
   standalone module

      A :term:`standalone module` contains functionality that is shared
      across scripts and other modules. Modules group functions and class
      definitions that address issues in resolving particular algorithmic
      problems.
   
Repository layout
=================

The repository contains the following directories:

:file:`CGATPipelines`
    ruffus_ pipeline scripts. Each pipeline consists of three parts.
    The pipeline itself is called :file:`pipeline_xyz.py`. The
    subdirectory :file:`pipeline_xyz` contains the default configuration
    files for the pipeline. The subdirectory
    :file:`pipeline_docs/pipeline_xyz` contains the CGATReport_ for the
    pipeline. Each pipeline is further associated with a file called
    :file:`PipelineXYZ.py` that contains utility function for that
    pipeline.

:file:`scripts`
    Utility scripts for managing and running pipelines.

:file:`doc`
    The sphinx_ documentation of the code repository.

:file:`tests`
    Testing code for pipeline modules.

:file:`R`
    R scripts and modules supporting pipeline functionality.


API Reference
=============

.. toctree::
   :maxdepth: 2

   scripts.rst
   modules.rst
   glossary.rst   

