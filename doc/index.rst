.. _contents:

==================================================================
CGATPipelines |version| - NGS and Genomics pipelines
==================================================================

This document brings together the various pipelines and scripts
written before and during CGAT. This documentation has two
parts. Below is the general documentation covering the complete
code collection. 

Overview
========

The CGAT pipeline collection has grown out of the work in comparative
genomics by the `Ponting
<http://www.dpag.ox.ac.uk/team/group-leaders/chris-ponting>`_ group in
the last decade. Now, CGAT_ has added functionality for
next-generation sequencing analysis (NGS).

CGAT pipelines perform basic tasks, are fairly generic and might be of
wider interest.

.. _cgatpipelines:

Active Pipelines
================

The pipelines in this section are being used and maintained.

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

Other Pipelines
===============

.. toctree::
   :maxdepth: 1

   CGATPipelines.rst

Background
==========

.. toctree::
   :maxdepth: 1

   InstallingPipelines.rst
   UsingPipelines.rst
   BuildingPipelines.rst
   PipelineReports.rst
   PipelinesBackground.rst   
   Release.rst

For information on how to contribute to the pipeline collection,
please see the `CGAT code collection
<https://www.cgat.org/downloads/public/cgat/documentation/>`_.

Reference
=========

This section contains a reference to the complete contents
of the pipeline collection.

.. toctree::
   :maxdepth: 1

   Reference.rst

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

