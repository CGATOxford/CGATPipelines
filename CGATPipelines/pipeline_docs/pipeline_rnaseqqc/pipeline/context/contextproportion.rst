.. _contextproportion:

======================
Transcript proportions
======================

This page summarises the proportions of reads that map to different genomic contexts
e.g:
  * protein-coding transcripts
  * small rna transcripts (rRNA, tRNA, snoRNA, Mt-RNA, snRNA etc)
  * repeat elements (LINEs, SINEs etc) 

Summary::
  * Aims of this analysis
  * What inputs/outputs
  * How the results were generated
  * What you should expect
  * Good example
  * Bad example
  * Links to other examples and the reasons that lie behind them

  This should be a text description of what to expect from the figures on this page.  What
  are the take-home messsages, how should the figures be interpreted, etc


Mapping summary
===============

The following table lists all genomic contexts that reads map to. 

.. report:: RnaseqqcReport.MappingContext
   :render: table
   :force:

   Number of alignments that align in various genomic contexts


Ribosomal and Repetitive RNA expression
=======================================

Ribosomal RNA is one of the most abundant transcripts in a cell
and dominates RNASeq samples until it is removed. The following
plots examine the number of alignments to ribosomal and repetitive 
RNA. Repetetive RNA annotation is taken from the UCSC repeatmasker 
tracks.

   Proportion of alignments that align to ribosomal annotations across samples.

.. report:: RnaseqqcReport.MappingContentRNA
   :render: stacked-bar-plot
   :split-at: 10


Context stats RNA expression
============================

The following plots examine the number of alignments to long intergenic 
non-coding and micro RNA. RNA annotations are taken from the 
UCSC repeatmasker tracks.

.. report:: RnaseqqcReport.MappingAltContent
   :render: table

   Number of alignments that align to long intergenic non-coding and 
   micro RNAs (from the UCSC repeatmasker track).

.. report:: RnaseqqcReport.MappingAltContent
   :render: stacked-bar-plot
   :split-at: 10

   Proportion of alignments that align to long intergenic non-coding and 
   micro RNAs across samples.


Aims
----

Inputs
------

Outputs
-------


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

   Add a comment about the bad example.  What is specifically bad about this example

More bad examples `<http://myBadData.html >`

Your data:



Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

