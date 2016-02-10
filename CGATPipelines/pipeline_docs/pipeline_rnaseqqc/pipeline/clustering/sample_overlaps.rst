.. _sample_overlaps:

=============================
Overlaps of sample expression
=============================

This page demonstrates the overlap between samples.  This is based on the intersection of
genes that are expressed with TPM > 1.  It can be used to identify obvious outliers,
and check that samples from the same conditions are more similar than between
conditions.

.. report:: RnaseqqcReport.SampleOverlapsExpress
   :render: sb-heatmap-plot
   :palette: GnBu

   Sample overlap of expressed genes with TPM > 1
