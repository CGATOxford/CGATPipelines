.. _memchip_details:

=============================
MEME-ChIP Detailed Results
=============================

The following is a list of the seed motifs found in each cluster:


.. report:: Motifs.MemeChipResults
   :render: table
   :groupby: track
   :large-html-class: sortable
   :force:

   Detailed MEME-ChIP results


.. report:: Motifs.CompFullSubsetLength
   :render: r-ggplot
   :tracks: r(_memechip)
   :groupby: track
   :layout: column-3
   :statement: aes(value) + geom_bar() + scale_x_log10() + facet_grid(slice~., scale="free_y") + theme_bw()

   Length distribution of full interval set vs subset used for meme motif calling
