.. _rnaseqqcpipeline:

===================
Top Expressed Genes
===================

Summary
=======
The top expressed genes should be similar for samples from the sample group.


Results
=======

.. report:: RnaseqqcReport.TopGenes
   :groupby: track
   :render: r-ggplot
   :statement: aes(gene_id, sample_name, fill=value) +
	       geom_tile() +
	       scale_fill_gradient2(low="goldenrod1", high="darkorchid4", mid="grey95", name="Z-score") +
	       xlab("") + ylab("") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(angle=90, vjust=0.5,
	       hjust=1)) +
	       facet_grid(factor_value~., scales="free")


Aims
----

Compare the expression of the top 20 most highly expressed genes across the samples

Inputs
------

Sailfish expression estimates (including bootstraps) from a single
sample subsampled to varying degrees (10-100% depth). Transcripts with
<1 average counts per sample are not included in this analysis.

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

