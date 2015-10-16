==========
CpG Tables
==========

Methylation values are provided in tabulated form below.
Samples are labelled as Tissue-Treatment-Replicate

Table columns are:
Tissue-Treatment-Replicate-perc = percentage methylation
Tissue-Treatment-Replicate-meth = number of methylated CpGs
Tissue-Treatment-Replicate-unmeth = number of unmethylated CpGs
contig, position, strand = genomic position of CpG
read_position = position of CpG on rrbs read
cpgi = Non-CpGIsland or CpGIsland
Tissue-Treatment-mean = mean methylation across replicates
Tissue-Treatment-var = methylation variance across replicates

Note: the read_position is calculated from an in silico MspI genome
digestion. A small fraction of reads will align outside of MspI
fragments. For these reads, the CpGs will not have their read_position
annotated

.. report:: Coverage.coveredCpGs
   :render: table
   :large: xls
   :force:
