.. _sample_overlaps:

=============================
Overlaps of sample expression
=============================

Summary
=======

The overlaps of sample expression provides a naive comparison between samples


Your Results
============


.. report:: RnaseqqcReport.imagesTracker
   :render: gallery-plot
   :glob: sailfish.dir/plots.dir/top_expressed_*png
	  

   Overlap between the top 1000 most highly
   expressed genes in each sample. Colour represents the fraction of
   genes in the intersection between the top 1000 genes in each pair
   of sample


About this section
==================

Aims
----
The overlaps of sample expression highlights the common overlapping
genes found between each of the samples. The aim of this page is to
identify any obvious outliers within your data and confirm that samples
from the same conditions are more similar than samples in different conditions.


Inputs
------
Sailfish expression estimates (TPM; transcripts per million) for each
transcriptstored in the 'transcript_quanitification' table within the pipeline csvdb database

Outputs
-------
A heatmap representing the overlapping (common) genes within the top
1000 genes in each sample

How the results are generated
-----------------------------

   1. For each sample, the top 1000 most highly expressed genes are
      identified
   2. An all-vs-all comparison of the samples is perfomed and the
      fraction of intersecting genes is rendered into a heatmap
      
   Note: Data manipulation and rendering are performed at the pipeline
   instead of reporting stage.

Examples
--------

   * Good example
   * Bad example
   * Links to other examples and the reasons that lie behind them

The good

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the good example.  What represents good data?

The bad

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad
   about this example

More bad examples `<http://myBadData.html >`



#Sample overlap - between samples
#--------------------------------

#.. report:: RnaseqqcReport.SampleOverlapsExpress
#   :render: sb-heatmap-plot
#   :title: Sample Overlap Heatmap

#   Sample Overlap Heatmap.
# NOTE: "nolabel-rows, nolabel-col" were not working!
#   :nolabel-cols: True
#   :nolabel-rows: True
