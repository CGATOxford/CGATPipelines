=====================================================================
Bias analysis results - Split by third identifier
=====================================================================

This page presents the analysis of potential biasing factors using
linear regression. The plot aesthetics are split by the third
identier, e.g replicate etc.

Genes/transcripts are binned according to their value for each
potential biasing factor (e.g GC content), with each bin containing an
equal number of genes/transcripts.  The mean expression for the
genes/transcripts for each sample is then calculated for each
bin. This mean expression is plotted below, along with a linear
regression or local (loess) regression for each sample.


Linear Regression
=================

.. report:: RnaseqqcReport.BiasFactors
   :render: r-ggplot
   :statement: aes(y=as.numeric(mean_expression), x=bin, colour=id_3)+
	       geom_point()+
	       stat_smooth(aes(group=sample,colour=id_3),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned factor levels. Linear regression.

Local Regression
=================

.. report:: RnaseqqcReport.BiasFactors
   :render: r-ggplot
   :statement: aes(y=as.numeric(mean_expression), x=bin, colour=id_3)+
	       geom_point()+
	       stat_smooth(aes(group=sample,colour=id_3),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned factor levels. Local regression.


