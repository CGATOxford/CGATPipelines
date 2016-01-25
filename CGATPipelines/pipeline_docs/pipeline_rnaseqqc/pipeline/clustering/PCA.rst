.. _PCA:

============================
Principal Component Analysis
============================

Principal component analysis (PCA) is a way of visualising datasets to identify 
the strong patterns that account for the majority of the variation between the 
datasets. 

Summary::
  * Aims of analysis
  * What inputs/outputs
  * How the results were generated
  * What you should expect
  * Good example
  * Bad example
  * Links to other examples and the reasons that lie behind them

  This should be a text description of what to expect from the figures on this page.  What
  are the take-home messsages, how should the figures be interpreted, etc

The good

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the good example.  What represents good data?

The bad

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad about this example

More bad examples `<http://myBadData.html >`

PCA Plot
========

.. report:: RnaseqqcReport.samplePCA
   :render: r-ggplot
   :statement: aes(x=PC1, y=PC2) +
	       geom_point() +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=20),
	       axis.text.y=element_text(size=20))

   Plot of First (PC1) and second (PC2) principal components from pricipal component
   analysis showing the latent variables that explain most of the variance in the dataset. 


   Graphs and tables
   Code snippets used to generate graphs and tables

Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

