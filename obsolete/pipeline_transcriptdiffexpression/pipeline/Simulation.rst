==================
Simulation Results
==================

The transcripts expression may be quantified here using one or more of
Kallisto, Sailfish or Salmon.  Like all RNA-Seq quantification
methods, these tools are expected to show poorer accuracy of
quantification for transcripts which have little unique sequence. To
identify and flag these transcripts, we perform a simulation here
where we generate a known number of transcripts (a ground truth),
represented as transcripts per million transcript (TPM) and compare
this to the estimated TPM. By repeating this multiple times (default
x30) we can calculate the correlation between the ground truth and the
estimated TPM across the simulations. Samples with poor correlation
are flagged in the results table as these could generate false
positives in the differential testing analysis. In addition, we can
compare the total estimated TPM across all simulations to the total
ground truth to identify transcripts where the TPM values are
systematically under or over-estimated. These transcripts are also
flagged.


Ground Truth TPM vs. Estimated TPM
==================================

These plots compare the sum of ground truth TPMs across all
simulations to the sum of estimated TPMs. N.B The sum of ground truths
is not the same for each transcript as the number of simulated copies
per transcript is randomly sampled [0-10] for each simulation. Where
the absolute fold-difference between the ground truth and estimate is
greater than 2-fold, the transcript is flagged.

.. report:: Simulations.simulationCorrelationsTpm
   :groupby: track
   :render: r-ggplot
   :statement: aes(tpm, est_tpm) +
	       geom_point(size=1, alpha=0.5,
	       aes(col=abs(as.numeric(as.character(log2diff_tpm)))<1)) +
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
	       name="Accuracy", labels=c("Poor", "Good")) +
	       xlab("Ground truth") +
	       ylab("Estimated TPM") +
	       geom_abline(col="grey40", size=0.5) +
	       geom_abline(intercept=0, slope=2, col="turquoise3",
                           size=0.5, linetype=2) +
	       geom_abline(intercept=0, slope=1, col="turquoise3",
                           size=0.5, linetype=2)

    Ground Truth TPM vs. Estimated TPM


Summary table for transcripts flagged by difference between ground
truth and estimated TPM


.. report:: Simulations.simulationCorrelationsSummaryFoldTpm
   :render: table
   :force:

   Summary table


In the following plot, the log2 fold difference between the total
ground truths and the total estimated TPM is show at different
levels of transcript "uniqueness". Absolute log2 fold differences
greater than 1 are plotted at -1 or 1 respectively and shown in
red. The "uniqueness" of a transcript is quantified by identifying the
fraction of all 31nt kmers which are unique to the transcript, i.e not
found in any other transcript.

.. report:: Simulations.simulationCorrelationsTpm
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
	           as.numeric(as.character(log2diff_tpm_thres))) +
	       geom_jitter(
	         alpha=0.3,
                 position=position_jitter(width=0.3, height=0),
		 aes(shape=abs(as.numeric(as.character(log2diff_tpm_thres)))==1,
		 size=abs(as.numeric(as.character(log2diff_tpm_thres)))==1,
		 colour=abs(as.numeric(as.character(log2diff_tpm_thres)))==1)) +
	       scale_size_manual(guide=FALSE, values=c(0.5, 1)) +
 	       scale_shape_discrete(guide=FALSE) + 
	       scale_colour_manual(guide=FALSE, values = c("grey15","red")) + 
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(axis.title.x=element_text(size=15),
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(breaks=as.factor(seq(0,2,0.2))) + 
	       xlab("Fraction unique kmers") +
	       ylab("Log2 Fold difference (Estimated/Ground Truth)")
	       
   Full plot	   


Ground Truth Counts vs. Estimated Counts
========================================


The plots below are identical except that the counts per transcript
are shown instead of the TPM

.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: r-ggplot
   :statement: aes(log(read_count,2), log(est_counts,2)) +
	       geom_point(size=1, alpha=0.5,
	       aes(col=abs(as.numeric(as.character(log2diff_counts)))<1)) +
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
	       name="Accuracy", labels=c("Poor", "Good")) +
	       xlab("Ground truth") +
	       ylab("Estimated Counts") +
	       geom_abline(col="grey40", size=0.5) +
	       geom_abline(intercept=0, slope=2, col="turquoise3",
                           size=0.5, linetype=2) +
	       geom_abline(intercept=0, slope=0.5, col="turquoise3",
                           size=0.5, linetype=2)

    Ground Truth Count vs. Estimated Count


Summary table for transcripts flagged by difference between ground
truth and estimated Count


.. report:: Simulations.simulationCorrelationsSummaryFoldCount
   :render: table
   :force:

   Summary table


.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
	           as.numeric(as.character(log2diff_counts_thres))) +
	       geom_jitter(
	         alpha=0.3,
                 position=position_jitter(width=0.3, height=0),
		 aes(shape=abs(as.numeric(as.character(log2diff_counts_thres)))==1,
		 size=abs(as.numeric(as.character(log2diff_counts_thres)))==1,
		 colour=abs(as.numeric(as.character(log2diff_counts_thres)))==1)) +
	       scale_size_manual(guide=FALSE, values=c(0.5, 1)) +
 	       scale_shape_discrete(guide=FALSE) + 
	       scale_colour_manual(guide=FALSE, values = c("grey15","red")) + 
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
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

.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: r-ggplot
   :statement: aes(fraction_unique,
	           as.numeric(as.character(log2diff_counts_thres))) +
	       geom_point()
	       theme_bw() +
	       theme(axis.title.x=element_text(size=15),
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       xlab("Fraction unique kmers") +
	       ylab("Log2 Fold difference (Estimated/Ground Truth)")
	       
   Full plot no binning		   


Correlation between ground truth TPM and estimated TPM
======================================================

These plots show the correlation between ground truth and estimated
TPM for each transcript against the "uniqueness" of the
transcript. The "uniqueness" of a transcript is quantified by
identifying the fraction of all 31nt kmers which are unique to the
transcript, i.e not found in any other transcript. Transcripts with
less than 3 % unique kmers are flagged.

    Correlation vs Fraction Unique Kmers

.. report:: Simulations.simulationCorrelationsTpm
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
	           as.numeric(as.character(tpm_cor))) +
	       geom_jitter(size=0.5, alpha=0.3, col="grey15",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(limits=seq(0,0.1,0.01)) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation (Estimated TPM vs ground truth)")
	       
   Zoomed plot	       


.. report:: Simulations.simulationCorrelationsTpm
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
	           as.numeric(as.character(tpm_cor))) +
	       geom_jitter(size=1, alpha=0.25, col="grey30",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(breaks=seq(0,1,0.1)) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation (Estimated TPM vs ground truth)")

    Full plot


Summary table for transcripts flagged by low fraction unique kmers

.. report:: Simulations.simulationCorrelationsSummaryKmers
   :render: table
   :force:

   Summary table


The plots below are identical except that the counts per transcript
are shown

    Correlation vs Fraction Unique Kmers

.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
                   as.numeric(as.character(counts_cor))) +
	       geom_jitter(size=0.5, alpha=0.3, col="grey15",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(limits=seq(0,0.1,0.01)) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation (Estimated Counts vs ground truth)")
	       
   Zoomed plot	       


.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: r-ggplot
   :statement: aes(as.factor(fraction_bin),
	           as.numeric(as.character(counts_cor))) +
	       geom_jitter(size=1, alpha=0.25, col="grey30",
                           position=position_jitter(width=0.3,
			   height=0)) +
	       stat_boxplot(geom ="errorbar", col="chartreuse3") +
	       geom_boxplot(outlier.colour=NA, col="chartreuse3") +
	       theme_bw() +
	       theme(
	       axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5),
	       axis.title.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       aspect.ratio=1) +
	       scale_x_discrete(breaks=seq(0,1,0.1)) + 
	       xlab("Fraction unique kmers") +
	       ylab("Correlation (Estimated Counts vs ground truth)")

    Full plot


These are the full tables of simulation results. The tpm and
read_count values are the sum of ground truths for the tpm and read
counts for the transcript across all the simulations. The
est_count/est_tpm values are the sum of estimated counts/tpm for the
transcript across all the simulations. The log2diff_count/log2diff_tpm
values are the log2-fold difference between the ground truth and
estimates. The count_cor/tpm_cor values are the pearson correlation
coefficient for the correlation between ground truth and estimated values

.. report:: Simulations.simulationCorrelationsTpm
   :groupby: track
   :render: xls-table

   Transcripts per million based correlations   
	    
.. report:: Simulations.simulationCorrelationsCount
   :groupby: track
   :render: xls-table

   Count-based correlations
