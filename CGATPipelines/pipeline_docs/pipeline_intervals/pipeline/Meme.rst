=====================
Ab-inito motif finding
=====================

Several Algorithms can be used to detect motifs, ab-inito and each can be run in several modes.

Thus far MEME and DREME are implimented in both descrimative and non-descriminative modes.

To increase the signal/noise ratio and reduce running time, subsets of sequence can be selected. 
Intervals are sorted by peak strength, score or shuffled randomly, and then either the top n
sequences or a fixed fraction of the sequences taken.

It is also possible to set a fixed width for intervals.

These settings can be different for each algorithm used.

The following criteria were used for each alogrithm:

.. report:: Motifs.SequenceSubset
   :render: table
   

The following is a summary of the results of each motif finding run, with links to the more detailed
results, including an analysis of the sequences used. 

:ref:`memedetails`

.. report:: Motifs.MemeRuns
   :render: table
   :groupby: slice
   :force:

   Meme results.


.. report:: Motifs.DremeRuns
   :render: table
   :groupby: slice
   :force:

   Dreme results




.. toctree::
   :glob:
   :hidden: 

   *Gallery.rst
   MemeGallery.rst
