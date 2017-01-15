========================================================================
Summary plots of the results for EdgeR differential expression analysis
========================================================================

This section presents some basic data exploration plots from the
transcript expression tables. 


Heatmaps and hierachical clustering alignment based
===================================================

Samples are expected to cluster by treatment group. The following
plots show the hierachical clustering of the samples based on the 200
most highly expressed transcripts. The gene expression values are
shown in the heatmap below the dendogram

The following plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*transcripts*heatmap.png

   clustering plots

The following plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*genes*heatmap.png

   clustering plots

Heatmaps and hierachical clustering pseudocounting
==================================================

Samples are expected to cluster by treatment group. The following
plots show the hierachical clustering of the samples based on the 200
most highly expressed transcripts. The gene expression values are
shown in the heatmap below the dendogram

Kallisto
--------

The following plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*transcripts*heatmap.png

   clustering plots

The following plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*genes*heatmap.png

   clustering plots

Salmon
--------

The following plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*transcripts*heatmap.png

   clustering plots

The following plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*genes*heatmap.png

   clustering plots

Prinicipal components analysis alignment based
==============================================
Principal components analysis is useful to see whether the major
sources of variation within the data seperate the samples as
expected. For example, do the sample seperate by their treatment
group. Only the first 4 PCs are shown here. 

PC1 and PC2
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*transcripts*pc1_pc2.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*genes*pc1_pc2.png

   PCA plots

PC3 and PC4
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*transcripts*pc3_pc4.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*genes*pc3_pc4.png

   PCA plots

PCA over 6 components
---------------------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*transcripts*variance.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*featurecounts*genes*variance.png

   PCA plots

Prinicipal components analysis pseudoalignment
==============================================
Principal components analysis is useful to see whether the major
sources of variation within the data seperate the samples as
expected. For example, do the sample seperate by their treatment
group. Only the first 4 PCs are shown here. 


Kallisto
========

PC1 and PC2
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*transcripts*pc1_pc2.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*genes*pc1_pc2.png

   PCA plots

PC3 and PC4
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*transcripts*pc3_pc4.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*genes*pc3_pc4.png

   PCA plots

PCA over 6 components
---------------------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*transcripts*variance.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*kallisto*genes*variance.png

   PCA plots

Salmon
========

PC1 and PC2
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*transcripts*pc1_pc2.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*genes*pc1_pc2.png

   PCA plots

PC3 and PC4
-----------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*transcripts*pc3_pc4.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*genes*pc3_pc4.png

   PCA plots

PCA over 6 components
---------------------

This plot shows the transcripts:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*transcripts*variance.png

   PCA plots

This plot shows the genes:

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/edger*salmon*genes*variance.png

   PCA plots


