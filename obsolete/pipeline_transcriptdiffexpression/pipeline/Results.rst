======================================
Differential expression results tables
======================================

DE Transcripts results tables
=============================

The following presents the results tables for the differential
expression (DE) testing for the transcripts. Only transcripts which
were identified as significantly differentially expressed by greater
than 2-fold are shown in the tables below. For the full results see
the table download links at the bottom of the page

The results table contains the following columns:

* gene_name = Ensembl gene name
* gene_id = Ensembl gene identification
* transcript_id = Ensembl transcript identification
* transcript_biotype = Ensembl transcript biotype
* expression = average expression of transcript (Transcripts per
  million; log2)
* fold = fold change in expression between conditions
* log2_fold = log2 transformed fold change
* p_value = p-value for null hypothesis that fold change is zero
* p_value_adj = p-value adjusted for multiple testing burden
  (Benjamini-Hochberg adjusted)
* significant = passes threshold for signficance = 1
* flagged = transcript is flagged as being difficult to quantify with
  simulated data. This may be becuase the transcript contains little
  sequence unique to itself (low_kmer), or else the estimated counts
  in the simulation may be > 2-fold different to the true counts
  (poor_accuracy). Either way, DE transcripts which are flagged are
  suspect.


These are the significant DE transcripts per test

.. report:: Results.SleuthResultsSig
   :render: table
   :groupby: track
   :large: xls
   :force:


   Results from DE testing


These are the full results tables

.. report:: Results.SleuthResults
   :render: xls-table
   :groupby: track
   :force:


   Results from DE testing


These are the results tables including the group means and the
expression values for each individual sample. All transcripts with
unadjusted p-values < 0.05 are included.

.. report:: Results.SummarisedResults
   :render: xls-table
   :groupby: track
   :force:

   Summarised results tables


DE Genes results tables
=======================

Deseq2 results

.. report:: Results.Deseq2ResultsSig
   :render: xls-table
   :groupby: track
   :force:

   Summarised results tables

