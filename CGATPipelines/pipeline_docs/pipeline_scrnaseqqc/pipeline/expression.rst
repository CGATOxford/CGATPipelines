.. _expression:

======================
Expression Concordance
======================

Here we present the concordance in expression of the spike-in transcripts for the two quantification methods, by alignment of short reads
followed by counting of reads intersecting transcriptome annotation, and lightweight alignment using Sailfish/Salmon..  It is assumed that
the spike-in transcripts are the ERCC92 RNAs, and thus are prefixed with ERCC-.  If this is not the case there will be no data displayed here.

.. report:: ScQCTrackers.CompareSpikeInExpression
   :render: ScQCTrackers.ExpressionFacetGrid
   
   A facet-grid comparing expression quantification of featureCounts and Sailfish on spiked-in transcripts
