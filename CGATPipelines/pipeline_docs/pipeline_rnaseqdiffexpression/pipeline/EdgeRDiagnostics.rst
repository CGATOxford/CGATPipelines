============================================================
Diagnostic plots for EdgeR differential expression analysis
============================================================

This section includes various diagnostic plots which are produced
during the differential transcript and gene expression analysis with EdgeR.

Results disprsion plots for alignment based counts
==================================================

The following plots the dispersion of the model that is fitted to the transcripts
following counting with featurecounts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*transcripts_results_dispersion.png

The following plots the dispersion of the model that is fitted to the genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*genes_results_dispersion.png

Results disprsion plots for pseudocounting
========================================

Kallisto
--------

The following plots the dispersion of the model that is fitted to the transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*transcripts_results_dispersion.png

The following plots the dispersion of the model that is fitted to the genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*genes_results_dispersion.png

Salmon
------

The following plots the dispersion of the model that is fitted to the transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*transcripts_results_dispersion.png

The following plots the dispersion of the model that is fitted to the genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*genes_results_dispersion.png


MA plots for alignment based counts
===================================
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for transcripts. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*transcripts_results*_MA_plot.png

   MA plots
	  
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for genes. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*genes_results*_MA_plot.png


   MA plots

MA plots for pseudocounting
===========================

Kallisto
--------
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for transcripts. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*transcripts_results*_MA_plot.png

   MA plots
	  
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for genes. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*genes_results*_MA_plot.png

   MA plots

Salmon
------
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for transcripts. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*transcripts_results*_MA_plot.png

   MA plots
	  
MA plots show the per-transcript mean abundance vs. the fold-change
between conditions for genes. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*genes_results*_MA_plot.png

   MA plots

Volcano plots for alignment based counts
========================================
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test for genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*transcripts*volcano_plot.png

   volcano plots

volcano plots show the per-gene fold change and p-value from the
differential expression statistical test for transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/featurecounts*gene*volcano_plot.png

   volcano plots


Volcano plots for pseudocounting
================================

Kallisto
--------
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test for genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*transcripts*volcano_plot.png

   volcano plots

volcano plots show the per-gene fold change and p-value from the
differential expression statistical test for transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/kallisto*gene*volcano_plot.png

   volcano plots

Salmon
------
volcano plots show the per-transcript fold change and p-value from the
differential expression statistical test for genes

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*transcripts*volcano_plot.png

   volcano plots

volcano plots show the per-gene fold change and p-value from the
differential expression statistical test for transcripts

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: DEresults.dir/edger/salmon*gene*volcano_plot.png

   volcano plots