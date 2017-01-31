=============================
Transcript expression summary
=============================

Expressed transcripts
=====================
Not all transcripts are expressed in all cells. Transcripts will also
be expressed at different levels. The following summarises the number
of transcripts expressed in each sample. Where the samples are
reasonably homogenous, we expect a similar number of expressed
transcripts per sample.

.. report:: Isoform.TranscriptExpressionOrdered
   :render: r-ggplot
   :statement: aes(y=index, x=value, colour=variable, group=variable) +
	       geom_line() +
	       theme_bw() +
	       theme(aspect.ratio=1,
	       axis.text.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.title=element_text(size=15),
	       legend.text=element_text(size=15)) +
	       xlab("Expression (Transcripts per Million (Log10))") +
	       ylab("Count") +
	       scale_colour_discrete(name="Group") +
	       ggtitle("Transcript expression")

   Distribution of expressed transcript per sample	

For most RNA-Seq experiments, we would expect to see the majority of
expressed transcripts being present in all samples. If a large number
of transcripts are only observed in a single sample, this may indicate
these transcripts are biological or technical noise. Alternatively,
the samples may be truly biologically heterogeneous. The plots below
shows the number of transcripts (y-axis) which are observed in the
number of samples (x-axis) at different expression thresholds
(Transcripts Per Million (log10))

.. report:: Isoform.TranscriptNumberSamplesExpressed
   :render: r-ggplot
   :statement: aes(x=as.factor(Count), y=No_transcripts) +
	       geom_line(aes(group=threshold, colour=as.factor(threshold))) +
	       theme_bw() +
	       theme(aspect.ratio=1,
	       axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.title=element_text(size=15),
	       legend.text=element_text(size=15)) +
	       xlab("# Samples") +
	       ylab("# Transcripts") +
	       scale_colour_discrete(name="Expression threshold") +
	       ggtitle("Transcript expression")

   Overlap of expressed transcripts	


Expression Correlation
======================
Depending on the heterogeneity of the samples, we might expect an
overall high correlation between transcript expression values. The
plot below presents all pairwise comparisons in the lower diagonal and
the absolute pearson correlation in the upper diagonal.

.. report:: Isoform.imagesTracker
   :render: gallery-plot
   :glob: summary_plots/*_pairwise_correlations.png
