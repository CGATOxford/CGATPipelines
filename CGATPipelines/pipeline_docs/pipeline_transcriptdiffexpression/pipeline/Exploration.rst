======================================
Transcript expression data exploration
======================================

This section presents some basic data exploration plots from the
transcript expression tables. 


Heatmaps and hierachical clustering
===================================
Samples are expected to cluster by treatment group. The following
plots show the hierachical clustering of the samples based on the 200
most highly expressed transcripts. The gene expression values are
shown in the heatmap below the dendogram

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/*heatmap.png

   clustering plots


Prinicipal components analysis
==============================
Principal components analysis is useful to see whether the major
sources of variation within the data seperate the samples as
expected. For example, do the sample seperate by their treatment
group. Only the first 4 PCs are shown here. 

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/*pc*_pc*.png

   PCA plots

