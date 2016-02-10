.. _sample_overlaps:

=============================
Overlaps of sample expression
=============================

This page demonstrates the overlap between samples.  This is based on the intersection of
genes that are expressed with TPM > 1.  It can be used to identify obvious outliers,
and check that samples from the same conditions are more similar than between
conditions.

This report presents a heatmap to show genes that are commonly expressed between samples.

Summary::
* Aims of this analysis: Identify which samples are expressing the same set of genes.
* What inputs/outputs: Input is transcription quantification per transcript and sample. This analysis uses TPM values with a threshold of TPM > 100. The output is a heatmap showing number of common genes between samples.
* How the results were generated: Extracted TPM information from SQLite database containing transcription quantification. Run a pairwise comparison to compute number of genes in common.

# NOTE: "nolabel-rows, nolabel-col" were not working!

.. report:: RnaseqqcReport.SampleOverlap
   :render: sb-heatmap-plot
   :title: Sample Overlap Heatmap
   :nolabel-cols: True
   :nolabel-rows: True

   Sample Overlap Heatmap.
