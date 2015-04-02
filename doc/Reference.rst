=========
Reference
=========

This section describes the layout of the code repository and contains
a reference to the complete contents of the pipeline collection.

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

Contents
========

.. toctree::
   :maxdepth: 1

   scripts.rst
   modules.rst
   glossary.rst   

