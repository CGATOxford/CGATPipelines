============================================================
Diagnostic plots for Sleuth differential expression analysis
============================================================

This section includes various diagnostic plots which are produced
during the differential transcript and gene expression analysis with Sleuth.

Results disprsion plots for pseudocounting
========================================

Kallisto
--------

The following plots the dispersion of the model that is fitted to the transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*transcripts_results_dispersion.png

The following plots the dispersion of the model that is fitted to the genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*genes_results_dispersion.png

Salmon
------

The following plots the dispersion of the model that is fitted to the transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*transcripts_results_dispersion.png

The following plots the dispersion of the model that is fitted to the genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*genes_results_dispersion.png


MA plots for pseudocounting
===========================

Kallisto
--------
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for transcripts. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*transcripts_results*_MA_plot.png

   MA plots
	  
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for genes. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*genes_results*_MA_plot.png

   MA plots

Salmon
------
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for transcripts. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*transcripts_results*_MA_plot.png

   MA plots
	  
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for genes. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*genes_results*_MA_plot.png

   MA plots


Volcano plots for pseudocounting
================================

Kallisto
--------
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test for genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*transcripts*volcano_plot.png

   volcano plots

volcano plots show the per-gene fold change and p-value from the
differential expression statistical test for transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/kallisto*gene*volcano_plot.png

   volcano plots

Salmon
------
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test for genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*transcripts*volcano_plot.png

   volcano plots

volcano plots show the per-gene fold change and p-value from the
differential expression statistical test for transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/sleuth/salmon*gene*volcano_plot.png

   volcano plots