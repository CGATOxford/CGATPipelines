===============
Genomic context
===============

This part of the report examines how tags are distributed
within the genome. The genome is divided into regions of
specific context, such as exon, intron, repeat, etc.

The counting is naive:

* Some genomic contexts can be overlapping, thus some tags will
  be counted several times.

* A tag needs to map to a context over at least 50% of its
  bases.  Thus some tags spanning several contexts might be
  dropped.

Global analysis
===============

All regions
===========

The following table lists all the genomic contexts that tags map to. 

.. report:: GenomicContext.ContextSummary
   :render: table
   :force:

   Number of tags that align in a certain genomic context

Significance
------------

The plots below show enrichments or depletions for particular
genomic annotations. Those which are signifant (FDR=5%) are shown
in red.

.. report:: GenomicContext.GatLogFoldContext
   :render: interleaved-bar-plot
   :groupby: track
   :colour: colour
   :layout: column-3
   :width: 300

   Log fold change. Those that are significant (qvalue < 0.05) are shown
   in red.

.. report:: GenomicContext.GatSignificantResultsContext
   :render: table
   :force:

   Significant gat results at an FDR of 5%
