.. _rnaseqqcpipeline:

===============================
Expression of genes of interest
===============================

Summary
=======
The top expressed genes should be similar for samples from the sample group.


Results
=======

.. report:: RnaseqqcReport.GenesOfInterest
   :groupby: track
   :render: r-ggplot
   :statement: aes(gene_id, sample_name, fill=log(TPM)) +
	       geom_tile() +
	       scale_fill_continuous(name="log TPM") +
	       theme_bw() +
	       xlab("") + ylab("") +
	       facet_grid(factor_value~., scales='free')


Aims
----

Compare the expression of particular genes of interest

Inputs
------

Sailfish expression estimates (Transcripts Per Million; TPM)

Outputs
-------

A heatmap of gene expression for each group of samples


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


Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

