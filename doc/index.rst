.. _contents:

==================================================================
CGATPipelines |version| - NGS and Genomics pipelines
==================================================================

This document brings together the various pipelines and scripts
written before and during the `CGAT Training Programme`_. This
documentation has two parts. Below is the general documentation
covering the complete code collection.

Overview
========

This is the code collecton written by members of the `CGAT Training Programme`_
for the execution of pipelines for next-generation sequencing analysis (NGS).

.. _cgatpipelines:

Documentation
=============

.. toctree::
   :maxdepth: 2

   PipelinesBackground.rst
   InstallingPipelines.rst
   UsingPipelines.rst
   Tutorials.rst
   BuildingPipelines.rst
   PipelineReports.rst
   Reference.rst
   Release.rst
   styleguide.rst
   Developers.rst

For information on how to contribute to the pipeline collection,
please see the `CGAT code collection
<https://www.cgat.org/downloads/public/cgat/documentation/>`_.

Active Pipelines
================

The pipelines in this section are being used and maintained by CGAT. We have travis
`integration <https://travis-ci.org>`_ testing of all of our scripts and we run
full local testing of our pipelines on `jenkins <https://jenkins.io>`_ each evening.

NGS Pipelines
-------------

.. toctree::
   :maxdepth: 1	

   pipeline_readqc: Read QC and processing before mapping <pipelines/pipeline_readqc.rst>
   pipeline_mapping: Short read mapping and QC <pipelines/pipeline_mapping.rst>
   pipeline_rnaseqtranscripts: RNA-seq transcript construction <pipelines/pipeline_rnaseqtranscripts.rst>
   pipeline_rnaseqdiffexpression: RNA-seq differential expression <pipelines/pipeline_rnaseqdiffexpression.rst>
   pipline_rnaseqlncrna: RNA-seq lncRNA analysis <pipelines/pipeline_rnaseqlncrna.rst>
   pipeline_peakcalling: ChIP-seq Peak calling <pipelines/pipeline_peakcalling.rst>
   pipeline_windows: Genomic read distribution <pipelines/pipeline_windows.rst>
   pipeline_exome: Exome-seq <pipelines/pipeline_exome.rst>

Genomics Pipelines
------------------

.. toctree::
   :maxdepth: 1	

   pipeline_annotations: Building genomic annotations <pipelines/pipeline_annotations.rst>
   pipeline_intervals: Annotating genomic intervals <pipelines/pipeline_intervals.rst>
   pipeline_ancestral_repeats: Deriving ancestral repeats <pipelines/pipeline_ancestral_repeats.rst>
   pipeline_chains: Processing genomic alignments <pipelines/pipeline_chains.rst>
   pipeline_liftover: Mapping genomic features <pipelines/pipeline_liftover.rst>

Meta pipelines
--------------

.. toctree::
   :maxdepth: 1	

   pipeline_testing: Regression testing of pipelines <pipelines/pipeline_testing.rst>

Disclaimer
==========

This collection of pipelines is the outcome of 10 years working in various 
fields in bioinformatics. It contains both the good, the bad and the ugly. 
Use at your own risk.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

