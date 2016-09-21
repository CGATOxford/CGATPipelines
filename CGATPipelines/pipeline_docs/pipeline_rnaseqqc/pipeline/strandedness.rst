.. _strandedness:

============
Strandedness
============

Summary
=======

Most RNA-Seq library preparations are now stranded. 


Results
=======

.. report:: RnaseqqcReport.PicardStrandBias
   :render: table
   :force:

.. report:: RnaseqqcReport.PicardStrandBias
   :render: r-ggplot
   :statement: aes(track, correct_strand_perc,
	       fill=correct_strand_perc>95) +
	       geom_bar(stat="identity") +
	       xlab("") +
	       ylab("% correct_strand") +
	       scale_fill_discrete(name=">95% correct strand") +	       
	       theme_bw() +
	       theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))


Aims
----

Check whether stranded library preparations are generating stranded
reads as expected. 


Inputs
------

Output from Picard RNA Seq QC tool


Outputs
-------

A table and bar plot with the percentage of reads which align to
annotated features with the correct strandedness. 

How the results are generated
-----------------------------

1. The CORRECT_STRAND_READS and INCORRECT_STRAND_READS columns are
   extracted from the picard_rna_metrics table and used to calculate
   the percentage of correctly stranded reads
2. Table and bar plots are rendered

What you should expect
----------------------

With a stranded RNA-Seq library preparation, the majority of reads
aligning to annotated features should do so on the same strand as the
feature. However, antisense reads will be generated from
almost all features to some extent, so one should not expect a 'correct
strand' value of 100%, even if the stranded library preparation worked
perfectly. A threshold of 95% is suggested for correct strandedness.


Examples
--------

  * Good example
  * Bad example
  * Links to other examples and the reasons that lie behind them

The good

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics


The bad

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics


More bad examples `<http://myBadData.html >`


Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.
