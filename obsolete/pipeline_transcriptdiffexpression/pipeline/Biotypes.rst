=======================================
Differential expression isoform biotype
=======================================

Isoform biotype
===============
Ensembl annotates all transcripts (isoforms) with a biotype (see:
http://www.ensembl.org/Help/Faq?id=468 ). Here, we examine the
differentially expressed isoforms for enriched biotypes.

The following plot presents each differential expression test in turn,
seperated by over- and under-expressed transcripts. The counts of each
biotype are shown within the plots:

Over-expressed transcripts:

.. report:: Isoform.TranscriptBiotypeSummaryUp
   :render: r-ggplot
   :statement: aes(x=significant, y=Count, fill=transcript_biotype) +
	       geom_bar(stat="identity", position="fill") +
	       geom_text(stat = "identity", lineheight=.8,
	       fontface="bold", aes(y = cumsum_centres,
	       label=paste0(Count, "\n(", round(100*fraction,1), "%)"))) +
	       theme_bw() +
	       theme(axis.text.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.title=element_text(size=15),
	       legend.text=element_text(size=15)) + xlab("") +
	       scale_fill_discrete(name="Transcript biotype") +
	       ylab("Fraction") +
	       ggtitle("Over-expressed")

   Over-expressed isoform biotype plots	       

Under-expressed transcripts:

.. report:: Isoform.TranscriptBiotypeSummaryDown
   :render: r-ggplot
   :statement: aes(x=significant, y=Count, fill=transcript_biotype) +
	       geom_bar(stat="identity", position="fill") +
	       geom_text(stat = "identity", lineheight=.8,
	       fontface", aes(y = cumsum_centres,
	       label=paste0(Count, "\n(", round(100*fraction,1), "%)"))) +
	       theme_bw() +
	       theme(axis.text.x=element_text(size=15),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.title=element_text(size=15),
	       legend.text=element_text(size=15)) + xlab("") +
	       scale_fill_discrete(name="Transcript biotype") +
	       ylab("Fraction") +
	       ggtitle("Under-expressed")

    Under-expressed isoform biotype plots	       


