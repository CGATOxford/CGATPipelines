==================
Inactive pipelines
==================

The pipelines are currently not being actively used. This might be
because they have evolved into different pipelines. For example,
pipeline_chipseq is now split into pipeline_mapping,
pipeline_peakcalling and pipeline_intervals. In other cases, pipelines
have addressed a specific issue but have not been reused since.

The pipelines are listed below for completeness:

Pipelines in development
========================

.. toctree::
   :maxdepth: 1	

   pipelines/pipeline_fusion.rst
   pipelines/pipeline_benchmark_rnaseqmappers.rst
   pipelines/pipeline_cufflinks_optimization.rst
   pipelines/pipeline_mappability.rst
   pipelines/pipeline_fastqToBigWig.rst
   pipelines/pipeline_mapping_benchmark.rst
   pipelines/pipeline_capseq.rst
   pipelines/pipeline_expression.rst
   pipelines/pipeline_transcriptome.rst
   pipelines/pipeline_polyphen.rst
   pipelines/pipeline_variant_annotation.rst
   pipelines/pipeline_variants.rst
   pipelines/pipeline_promotors.rst
   pipelines/pipeline_transfacmatch.rst
   pipelines/pipeline_exome_cancer.rst
   pipelines/pipeline_genesets.rst
   pipelines/pipeline_idr.rst
   pipelines/pipeline_metagenomecommunities.rst
   pipelines/pipeline_motifs.rst
   pipelines/pipeline_proj020.rst
   pipelines/pipeline_rnaseqqc.rst
   pipelines/pipeline_rrbs.rst
   pipelines/pipeline_timeseries.rst


Obsolete pipelines
==================

.. toctree::
   :maxdepth: 1	

   pipelines/pipeline_rnaseq.rst
   pipelines/pipeline_chipseq.rst
   pipelines/pipeline_medip.rst

=================
Lecgacy pipelines
=================

Within of the Ponting group we have developed a few other pipelines
that are based on the GNU ``make`` utility:

.. toctree::
   :maxdepth: 1

   make_pipelines/CompareTranscripts.rst
   make_pipelines/Gpipe.rst
   make_pipelines/MapTranscripts454.rst
   make_pipelines/Optic.rst

Some of these pipelines are still being used, though they are not actively
supported any more.
