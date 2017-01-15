.. _saturation:

==================================
Sequencing Depth Saturation curves
==================================

Summary
=======

If all other variables remain constant, increased sequencing depth
will increase quantification accurary. However the gains will be
minimal where sequencing depth is already sufficient. Here, we inspect
the impact of reducing the sequencing depth on the quantification
accuracy, using the expression estimates from the full depth sample as
a proxy "ground truth". If we can reduce the read depth to 50% without
having a serious impact on the quantification estimates, this suggests
the sequencing depth is sufficient.


Results
=======

Sequencing depth vs. quantification accuracy
============================================

.. report:: RnaseqqcReport.imagesTracker
   :render: gallery-plot
   :glob: sailfish.dir/plots.dir/saturation_plots_accuracy.png

   Average absolute difference between the expression estimate and the
   expected expression, normalised to the expected expression. The
   lines represent deciles of gene expression


Sequencing depth vs. boostrap variance
======================================

.. report:: RnaseqqcReport.imagesTracker
   :render: gallery-plot
   :glob: sailfish.dir/plots.dir/saturation_plots_boostrap_cv.png

   Average CV per sequencing depth for deciles of expression since
   expression level

Aims
----

Assess whether the sequencing depth is sufficient

Inputs
------

Sailfish expression estimates (including bootstraps) from a single
sample subsampled to varying degrees (10-100% depth). Transcripts with
<1 average counts per sample are not included in this analysis.

Outputs
-------

Two figures are generated to assess the impact of sequencing depth on
quantification accruacy.

1. Sequencing depth vs. quantification accuracy

Here we want to see how reducing the sequencing depth affects the
accuracy of quantification. We don't have a ground truth for the
expression estimate so the expression estimate from the full depth
file is used as the expected expression. The figure shows the average
absolute difference between the expression estimate and the expected
expression, normalised to the expected expression. E.g if expected is
10 and observed is 6 the normalised absolute difference is:
|(6 - 10) / 10| = 0.4. The lines represent deciles of expression since
expression level is correlated with quantification accuracy

2. Sequencing depth vs. boostrap variance

Here we want to see how reducing the sequencing depth affects the
sailfish boostrap estimates. The variance in the boostrap estimates is
a proxy for the technical variance and should decrease with increased
sequencing depth. The boostrap variance is presented as the
coefficient of variance (CV; sd/mean) to make it comparable over
different expression levels.The lines present average CV per
sequencing depth for deciles of expression since expression level is
correlated with quantification accuracy


The good (both plots)

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   We would expect that the top expressed transcripts (e.g deciles 9 &
   10 ) should show very little difference (<0.1) in expression
   between 50-100% sequencing depth.

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   We would expect that the top expressed transcipts (e.g deciles 9 &
   10 ) should show very little difference (<0.1) in expression between
   50-100% sequencing depth.


The bad (both plots)

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad about this example

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad about this example


More bad examples `<http://myBadData.html >`






Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

