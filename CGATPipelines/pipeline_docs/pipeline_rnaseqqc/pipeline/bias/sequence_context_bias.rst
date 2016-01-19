.. _sequence_context:

=======================
Sequence context biases
=======================

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
genes/transcripts for each sample is then calculated for each
bin.

This mean expression is plotted below, along with a local (loess)
regression for each sample. The expectation is that the fit for each
individual sample will be very similar

The plot aesthetics are split in turn by the first, second and third
identier, e.g genotype, condition, replicate etc.

Following this, the relationship between the potential biasing factor
and expression level is explored by computing the Spearman rank
correlation. The gradient of the linear regression and rho value of
the correlation are shown. The expectation is that the gradient and
rho for each individual sample will be similar. Note: some of the
relationships are not linear. In these cases, the gradient may not be
meaningfull

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


Local Regression
=================

.. report:: RnaseqqcReport.BiasFactors
   :render: r-ggplot
   :statement: aes(y=as.numeric(mean_expression), x=bin, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=sample,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned factor levels. Local regression.

.. report:: RnaseqqcReport.BiasFactors
   :render: r-ggplot
   :statement: aes(y=as.numeric(mean_expression), x=bin, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=sample,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned factor levels. Local regression.

.. report:: RnaseqqcReport.BiasFactors
   :render: r-ggplot
   :statement: aes(y=as.numeric(mean_expression), x=bin, colour=id_3)+
	       geom_point()+
	       stat_smooth(aes(group=sample,colour=id_3),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Third Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned factor levels. Local regression.

   Graphs and tables
   Code snippets used to generate graphs and tables


Spearmons rho summary plots
===========================

.. report:: RnaseqqcReport.CorrelationSummaryGC
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryA
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryT
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryC
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryG
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.
    

Linear regression gradient summary plots
========================================

.. report:: RnaseqqcReport.GradientSummaryGC
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryA
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.


.. report:: RnaseqqcReport.GradientSummaryT
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryC
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryG
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.


Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

