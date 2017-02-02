.. _rand_dreme_details:

=================================
Randomized DREME Detailed Results
=================================

The following table lists all the motifs found.

.. report:: Motifs.DremeResults
   :render: table
   :groupby: track
   :slices: rand_dreme
   :force:

   Meme results, showing for each motif the width,
   the number of times it is found, its evalue
   and its information content.

The following is a comparision of the subset of sequences used for Motif finding compared to the full
interval set:


.. report:: Motifs.CompFullSubsetLength
   :render: r-ggplot
   :tracks: r(_dreme)
   :groupby: track
   :layout: column-3
   :statement: aes(value) + geom_bar() + scale_x_log10() + facet_grid(slice~., scale="free_y") + theme_bw()

   Length distribution of full interval set vs subset used for dreme motif calling

