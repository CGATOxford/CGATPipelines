======================================
Differential expression results tables
======================================

DE results tables
=================

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


These are the tables of Counts and TPM per gene per sample, with group means

#.. report:: Results.GeneLevelExpression
#   :render: xls-table
#   :groupby: track
#   :force:

   Transcripts per million per sample, aggregated per gene

These are the tables of Counts and TPM per transcript per sample, with
gene ID annotations and simulation flags.

.. report:: Results.SleuthTpmAll
   :render: xls-table
   :groupby: track
   :force:


   Transcripts per million per sample

.. report:: Results.SleuthCountsAll
   :render: xls-table
   :groupby: track
   :force:

   Counts per sample
