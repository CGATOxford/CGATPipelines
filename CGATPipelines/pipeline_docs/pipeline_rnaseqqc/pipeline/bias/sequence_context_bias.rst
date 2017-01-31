.. _sequence_context:

=======================
Sequence context biases
=======================

Summary
=======

The following presents the analysis of potential biasing
factors. Library preparation and sequencing technologies are known to
show biases against particular sequence contexts. For example,
illumina sequencing is biased against regions of extreme GC content. The
aim of this analysis is to establish whether the biases are consistent
between samples and particularly between groups of samples (e.g
different conditions)

Genes/transcripts were binned according to their value for each
potential biasing factor (e.g GC content), with each bin containing an
equal number of genes/transcripts.  The mean expression for the
genes/transcripts within each bin is calculated for each sample. This
mean expression is plotted below, along with a local (loess)
regression for each sample. The expectation is that the fit for each
individual sample will be very similar.


Results
=======

.. report:: RnaseqqcReport.BiasFactorsGCContent
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("GC content") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned GC content. Local regression.

.. report:: RnaseqqcReport.BiasFactorsAA
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("AA dinculeotides") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned AA dinucleotide content. Local regression.

.. report:: RnaseqqcReport.BiasFactorsCC
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("TT dinculeotides") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned CC dinucleotide content. Local regression.

.. report:: RnaseqqcReport.BiasFactorsGG
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("CC dinculeotides") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned GG dinucleotide content. Local
   regression.

.. report:: RnaseqqcReport.BiasFactorsTT
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("GG dinculeotides") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned TT dinucleotide content. Local
   regression.


.. report:: RnaseqqcReport.BiasFactorsLength
   :render: r-ggplot
   :statement: aes(bin, value, colour=factor_value)+
	       geom_point()+
	       stat_smooth(aes(group=sample_id, colour=factor_value),
	                       method=loess, se=F)+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       ggtitle("Length") +
	       theme_bw() +
	       theme(
	       aspect.ratio=1,
	       axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned length. Local regression.



Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.


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

Your data:

.. report:: RnaseqqcReport.BiasFactors
   :render: table


About this section
==================

Inputs
------
A table of expression values in TPM, binnned by various
possible biasing facors such as GC content.

Outputs
-------
Scatter plots showing the relationship between the potential biasing
factor and expression 


What you should expect from results
-----------------------------------

Sequencing biases will undoubtedly be present. For example, illumina
sequencing is known to be bias against extremes of GC content (< 20 %
/ > 80 %), so extreme GC content transcripts will appear to be more
lowly expressed. So long as the downstream analysis is relative, e.g
change in expression between two conditions, this need not be a
problem so long as the bias is consistent. However, if the bias is
inconsistent and the difference in the bias is confounded with the
experimental design, e.g condition 1 shows a distinct sequencing bias
from condition 2, then differential expression testing between the two
conditions may be invalid. To resolve this, you will need to consider
tools to adjust for the sequencing bias.
