.. _sample_overlaps:

=============================
Overlaps of sample expression
=============================

This report presents a heatmap to show genes that are commonly expressed between samples.

Summary::
  * Aims of this analysis: Identify which samples are expressing the same set of genes.
  * What inputs/outputs: Input is transcription quantification per transcript and sample. This analysis uses TPM values with a threshold of TPM > 100. The output is a heatmap showing number of common genes between samples.
  * How the results were generated: Extracted TPM information from SQLite database containing transcription quantification. Run a pairwise comparison to compute number of genes in common.
  * What you should expect: To-Do
  * Good example: To-Do
  * Bad example: To-Do
  * Links to other examples and the reasons that lie behind them: To-Do

  This should be a text description of what to expect from the figures on this page.  What
  are the take-home messsages, how should the figures be interpreted, etc

TO-DO: The good

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the good example.  What represents good data?

TO-DO: The bad

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad about this example

TO-DO: More bad examples `<http://myBadData.html >`

Your data:

# NOTE: "nolabel-rows, nolabel-col" were not working!

.. report:: RnaseqqcReport.SampleOverlap
   :render: sb-heatmap-plot
   :title: Sample Overlap Heatmap
   :nolabel-cols:
   :nolabel-rows:

   Sample Overlap Heatmap.

Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

