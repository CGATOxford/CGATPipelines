======================================
Differential expression results tables
======================================

DE results tables
=================

The following presents the results tables for the differential
expression (DE) testing for the transcripts. 

The results table contains the following columns:

* gene_name = Ensembl gene name
* gene_id = Ensembl gene identification
* transcript_id = Ensembl transcript identification
* transcript_biotype = Ensembl transcript biotype
* expression = average expression of transcript
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


These are the significant DE transcript

.. report:: Results.SleuthResultsSig
   :render: table
   :large: xls
   :force:


   Results from DE testing


This is the full results table

.. report:: Results.SleuthResults
   :render: xls-table
   :force:


   Results from DE testing
