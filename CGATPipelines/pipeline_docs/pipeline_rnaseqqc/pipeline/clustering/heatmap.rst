.. _heatmap:

=======
Heatmap
=======

This page illustrates the similarity between samples based on expression of all genes
in each sample using the Pearson correlation coefficient.  This can be used to identify
obvious outliers.  Samples are hierarchically clustered using average linkage
clustering.

.. report:: RnaseqqcReport.SampleHeatmap
   :render: sb-heatmap-plot
   :palette: GnBu

   Similarity of samples using pair-wise Pearson correlations
