============================================================
Diagnostic plots for Sleuth differential expression analysis
============================================================

This section includes various diagnostic plots which are produced
during the differential transcript expression analysis with Sleuth.

Techincal and Biological variance
=================================
Sleuth uses the bootstapping values from Kallisto to estimate the
per-transcript technical variance. The plots below show the relationship
between the estimates for the techincal vairance and total variance
(techincal + biological).


.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/*vars.png

   variance plots


Mean-Variance plots
===================
Sleuth expects a relationship to be present between the mean
expression of a gene and the variance in expression across the
dataset. These plots present the relationship. Each dot is a
transcript.

To more accurately estimate the true per-transcript variance, information is
shared between genes expressed at a similar level. The blue dots
represent transcripts which are used to model the relationship between
mean expression and variance. The red curve represents the smoothed
fit.

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/*mean_var.png

   mean-variance plots


MA plots
========
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/*MA_plot.png

   MA plots
	  

Volcano plots
=============
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/*volcano_plot.png

   volcano plots
