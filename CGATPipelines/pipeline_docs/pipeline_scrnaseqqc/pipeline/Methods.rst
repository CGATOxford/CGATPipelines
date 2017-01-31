=======
Methods
=======

A key step in any scientific experiment is to quality control the data that are generated.  This is to prevent spurious signals downstream
and increase the signal:noise ratio in the data.  This is especially important in single cell experiments that consider cell-to-cell
heterogeneity.  Technical noise is the overriding source of variation in single-cell experiments, thus reducing this by systematic and
extensive QC is key to maximising genuine biological signals in downstream analyses.

This pipeline collects QC measures derived from standard NGS workflows, such as read QC and read alignment statistics, as well as those
based on expression of native genes and synthetic spiked-in transcripts.

Following some optimisation and assessment of QC measures it was found that, generally speaking, poorly performing cells do not fall into
their own clusters.  Rather they are better represented by groups of outliers dependning on the QC measure that identifies them as being
of poor quality.  This is largely in line with the findings of Ilicic et al :pmid:`26887813`, where a number of different features are required
to identify poorly performing cells.  Their method is designed for cells captured using the Fluidigm C1 autoprep system, which does not apply
to cells FACS sorted directly into plates for library prep.

There are a number of reasons a cell may perform poorly, these are the ones I can think of::
  * dead cell
  * pre-apoptotic cell
  * empty well
  * cell lyses during sorting process
  * low RNA content of cell
  * poor/inefficient cDNA conversion during reverse transcription
  * sub-optimal reagent concentrations in well for library prep

The QC pipeline does not aim to identify what the explanation of poor quality cells, only to identify them to be excluded for downstream analysis.
This report goes over the marginal and joint distributions of the QC measures. It then performs a PCA to find the major axis of variation, on the assumption
that the greatest differences will be between good and poor quality cells.  A number of heurtistic thresholds are then combined to identify poorly
performing cells.  Another PCA is used to identify which QC variables best discriminate the good and poor performing cells.  These are also visualised on a
low dimensional representation of the cells through both the spiked-in transcripts and the detected protein-coding transcripts.

The report will highlight which cells specifically should be removed prior to downstream analysis, with a short description of what measure they fall
short on.

The measures considered are::
  * % of perfectly mapped pairs
  * nucleotide mismatch rate
  * strand balance
  * % duplicated reads
  * % reads overlapping to rRNA annotations
  * % reads overlapping to protein-coding transcripts
  * median insert size
  * median absolute deviation of insert size
  * % total expression from spike-in transcripts
  * number of genes expressed, TPM > 0

The report also compared the expression quantification output from the lightweight-alignment method (currently Sailfish) and read alignment based on
feature overlap (currently featureCounts).  It only considers the concordance in expression for the spike-in genes, and can demonstrate which are and
are not detected concordantly.  Discordance in the expression measurements may also be an indicator of poorly performing cells.

Data Visualisation
------------------
QC measure data are visualised with a combination of scatter, kernal density and bar plots.  These are designed to presnet what I perceive
as being the most important features of these QC data.  As with any high-throughput experimental data, metadata is also vital for the proper
interpretation of results.  It is expected that a table of metadata exists in the pipeline CSVDB, and that certain fields are present, namely
`Plate`, `seqRun` and `Well`.  Other columns may be present, but for now are ignored.  These metadata features are used to visualise the data
only, and it is down to the user to interpret their patterns and meaning.  The principal reason for their use is to aid in the identification
of obvious explanatory factors for poor quality data.  For instance there may be several plates of cells in a single sequencing run, one or more
of which may contain the lion's share of poor quality cells.  There may have been a pipetting of sorting problem with a section of the plate, which
would manifest in the wells, or a problem with one of the sequencing runs, etc. 
