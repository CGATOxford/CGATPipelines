=====================
Picard RNASeq Metrics
=====================

This section looks at the output from the Picard RNASeqMetrics.
If in doubt please check for the presence of a histogram in the ".....picard_rna_metrics" file(s) in the bam directory. 

Picard CollectRnaSeq metrics for each BAM file for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#CollectRnaSeqMetrics>`_
for a definition of the field contents.

Data table

.. report:: Mapping.PicardRnaMetrics
   :render: table

All histograms. In the plot the X axis the relative position in the transcript, where 0 is the 5' end and 100 represents the 3' end.
The y axis represents the relative coverage of that part of a transcript.

.. report:: Mapping.RnaBiasTable
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :xrange: 0,100
   :ytitle: read coverage
   :xtitle: 5' position 3'

   Individual plots

.. report:: Mapping.RnaBiasTable
   :render: line-plot
   :as-lines:
   :xrange: 0,100
   :yrange: 0,
   :groupby: track
   :ytitle: read coverage
   :xtitle: 5' position 3'
   :layout: column-3
   :width: 200

   Grouped plots

