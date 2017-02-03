.. _disc_meme_details:

====================================
Descriminative Meme Detailed Results
====================================

The following table lists all the motifs found.

.. report:: Motifs.MemeResults
   :render: table
   :groupby: track
   :slices: disc_meme
   :force:

   Meme results, showing for each motif the width,
   the number of times it is found, its evalue
   and its information content.

The following is a comparision of the subset of sequences used for Motif finding compared to the full
interval set:


.. report:: Motifs.CompFullSubsetLength
   :render: r-ggplot
   :tracks: r(_meme)
   :groupby: track
   :layout: column-3
   :statement: aes(value) + geom_bar() + scale_x_log10() + facet_grid(slice~., scale="free_y") + theme_bw()

   Length distribution of full interval set vs subset used for meme motif calling


The following is a comparison of the positive and negative sets used for the motif finding

.. report:: Motifs.PosNegMemeLength
   :render: box-plot
   :groupby: track
   :logscale: y
   :ytitle: length(bp)
   :layout: column-3

   Length distribution of positive and negative sequence sets used for peak calling


.. report:: Motifs.PosNegMemeGC
   :render: box-plot
   :groupby: track
   :ytitle: Proportion GC
   :layout: column-3

   Distribution of GC contents for positive and negatives sequence sets used for peak calling


.. report:: Motifs.PosNegMemeN
   :render: line-plot
   :transform: histogram
   :groupby: track
   :xtitle: Proportion of bases N
   :as-lines:
   :ytitle: Cumulative frequence
   :tf-aggregate: cumulative
   :tf-bins: 20

   Fraction of bases N (possibly masked) for positive and negative sequence sets used for peak calling

