==================
Simulation Results
==================

The transcripts expression is quantified here using a tool called
Kallisto.  Like all RNA-Seq quantification methods, Kallisto is
expected to show poorer accuracy of quantification for transcripts
which have little unique sequence. To identify and flag these
transcripts, we perform a simulation here where we generate a known
number of reads per transcript (a ground truth) and compare this to
the estimated counts from Kallisto. By repeating this multiple times
(default x30) we can calculate the correlation between the ground
truth and the estimated counts across the simulations. Samples with
poor correlation are flagged in the results table as these could
generate false positives in the differential testing analysis. In
addition, we can compare the total estimated counts across all
simulations to the total ground truth to identify transcripts where
the counts are systematically under or over-estimated. These
transcripts are also flagged.


Ground Truth vs. Estimated Counts
=================================

These plots compare the sum of ground truths across all simulations to
the sum of estimated counts. N.B The sum of ground truths is not the
same for each transcript as the number of simulated reads per
transcript is randomly sampled [0-100] for each simulation. Where the
absolute fold-difference between the ground truth and estimate is
greater than 1.5-fold, the transcript is flagged.

.. report:: Simulations.simulationCorrelations
   :render: r-ggplot
   :statement: aes(read_count, est_counts) +
	       geom_point(size=1, alpha=0.5,
	       aes(col=abs(as.numeric(as.character(log2diff)))<0.585)) +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.x=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.title=element_text(size=15),
	       legend.text=element_text(size=15),
	       aspect.ratio=1) +
	       scale_colour_discrete(
	       name="Kallisto\naccuracy", labels=c("Poor", "Good")) +
	       xlab("Ground truth") +
	       ylab("Estimated counts") +
	       geom_abline(col="grey40", size=0.5) +
	       geom_abline(intercept=0, slope=1.5, col="turquoise3",
                           size=0.5, linetype=2) +
	       geom_abline(intercept=0, slope=0.66, col="turquoise3",
                           size=0.5, linetype=2)

    Ground Truth vs. Estimated Counts


Summary table for transcripts flagged by difference between ground
truth and estimated counts

.. report:: Simulations.simulationCorrelationsSummaryFold
   :render: table
   :force:

   Summary table


In the following plot, the log2 fold difference between the total
ground truths and the total estimated counts is show at different
levels of transcript "uniqueness". Absolute log2 fold differences
greater than 1 are plotted at -1 or 1 respectively and shown in
red. The "uniqueness" of a transcript is quantified by identifying the
fraction of all 31nt kmers which are unique to the transcript, i.e not
found in any other transcript.

.. report:: Simulations.simulationCorrelations
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin), as.numeric(as.character(log2diff_thres))) +
	       geom_jitter(alpha=0.3,
                           position=position_jitter(width=0.3, height=0),
			   aes(shape=abs(as.numeric(as.character(log2diff_thres)))==1,
			   size=abs(as.numeric(as.character(log2diff_thres)))==1,
			   colour=abs(as.numeric(as.character(log2diff_thres)))==1)) +
	       scale_size_manual(guide=FALSE, values=c(0.5, 1)) +
 	       scale_shape_discrete(guide=FALSE) + 
	       scale_colour_manual(guide=FALSE, values = c("grey15","red")) + 
	       stat_boxplot(geom ='errorbar', col="chartreuse3") +
	       geom_boxplot(outlier.col=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(axis.title.x=element_text(size=15),
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(breaks=as.factor(seq(0,1,0.1))) + 
	       xlab("Fraction unique kmers") +
	       ylab("Log2 Fold difference (Estimated/Ground Truth)")
	       
   Full plot	   


Correlation between ground truth and estimated counts
=====================================================

These plots show the correlation between ground truth and estimated
counts for each transcript against the "uniqueness" of the
transcript. The "uniqueness" of a transcript is quantified by
identifying the fraction of all 31nt kmers which are unique to the
transcript, i.e not found in any other transcript. Transcripts with
less than 3 % unique kmers are flagged.

    Correlation vs Fraction Unique Kmers

.. report:: Simulations.simulationCorrelations
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin), as.numeric(as.character(cor))) +
	       geom_jitter(size=0.5, alpha=0.3, col="grey15",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ='errorbar', col="chartreuse3") +
	       geom_boxplot(outlier.col=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(limits=as.factor(seq(0,0.1,0.01))) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation")
	       
   Zoomed plot	       


.. report:: Simulations.simulationCorrelations
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin), as.numeric(as.character(cor))) +
	       geom_jitter(size=1, alpha=0.25, col="grey30",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ='errorbar', col="chartreuse3") +
	       geom_boxplot(outlier.col=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(breaks=as.factor(seq(0,1,0.1))) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation")

    Full plot


Summary table for transcripts flagged by low fraction unique kmers

.. report:: Simulations.simulationCorrelationsSummaryKmers
   :render: table
   :force:

   Summary table
