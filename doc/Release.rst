=============
Release notes
=============

The pipelines are under continuous improvement and the
latest code can always be found in the code repository.
Nevertheless, we occasionally prepare releases. Notes on
each release are below.

0.3.3
-----

* Updated documentation for pipeline_readqc tutorial; https://github.com/CGATOxford/CGATPipelines/pull/386
* Updated multiqc to 1.4 in conda environment ; https://github.com/CGATOxford/CGATPipelines/pull/389
* pipeline_mapping: merge bam files as part of the full target; https://github.com/CGATOxford/CGATPipelines/pull/384
* fix Cluster.py for PBSpro ; https://github.com/CGATOxford/CGATPipelines/pull/385
* pipeline_peakcalling: produce normalised BigWigs if the samples are generated using the quantitative ChIP-seq method ; https://github.com/CGATOxford/CGATPipelines/pull/382
* pipeline_readqc: bugfix trim_galore pre-processing ; https://github.com/CGATOxford/CGATPipelines/pull/391
* Control module in Pipeline.py: bugfix the PBSpro connection ; https://github.com/CGATOxford/CGATPipelines/pull/390
* pipeline_enrichment: bugfix runGsea task in pipeline_enrichment ; https://github.com/CGATOxford/CGATPipelines/pull/392
* Cluster module in Pipeline.py: bugfix the job submission in SGE with an empty queue parameter ; https://github.com/CGATOxford/CGATPipelines/pull/393
* pipeline_bamstats: changed way bam files are imported ; https://github.com/CGATOxford/CGATPipelines/pull/395
* pipeline_exome: fix import of PipelineExomeAncestry ; https://github.com/CGATOxford/CGATPipelines/pull/397
* pipeline_rnaseqdiffexpression: small change to avoid confusion when running featurecounts ; https://github.com/CGATOxford/CGATPipelines/pull/398/files
* pipeline.py: prevent pipelines from running when DRMAA is not available ; https://github.com/CGATOxford/CGATPipelines/pull/399
* make threads=1 instead of threads=10 for build_report in pipeline.ini ; https://github.com/CGATOxford/CGATPipelines/pull/401
* fix memory requirements for pipeline_genesets ; https://github.com/CGATOxford/CGATPipelines/pull/402
* update docs to clarify expected fastq name for input files ; https://github.com/CGATOxford/CGATPipelines/pull/403
* replace pipeline_annotations with pipeline_genesets in quickstart ; https://github.com/CGATOxford/CGATPipelines/pull/405
* update report for bamstats ; https://github.com/CGATOxford/CGATPipelines/pull/404
* add back explicit reporting of failed cluster jobs ; https://github.com/CGATOxford/CGATPipelines/pull/406
* use cgat gtf2tsv -f in pipeline_genesets ; https://github.com/CGATOxford/CGATPipelines/pull/408
* make sure we stick to conda 4.3 until a workaround is found for conda 4.4 ; https://github.com/CGATOxford/CGATPipelines/pull/409
* increase priority of $HOME/.cgat ini file ; https://github.com/CGATOxford/CGATPipelines/pull/412 ; https://github.com/CGATOxford/CGATPipelines/pull/414
* typo in pipeline_bamstats ; https://github.com/CGATOxford/CGATPipelines/pull/413
* make sure the job environment is the same when submitting jobs to the cluster or running them locally; https://github.com/CGATOxford/CGATPipelines/pull/411
* remove duplicated function ; https://github.com/CGATOxford/CGATPipelines/pull/416
* pipeline_exome: update snpEff from 4.1 to 4.3 ; https://github.com/CGATOxford/CGATPipelines/pull/417
* separate statements for paired-end and single-end data in CGATPipelines/PipelineRnaseq.py ; https://github.com/CGATOxford/CGATPipelines/pull/418
* pipeline_mapping: avoid duplicating entries in loadReadCounts ; https://github.com/CGATOxford/CGATPipelines/pull/420
* added pipeline_chiptools ; https://github.com/CGATOxford/CGATPipelines/pull/379
* fix cutadapt paired untrimmed command ; https://github.com/CGATOxford/CGATPipelines/pull/423
* pipeline_readqc: disable fastq_screen by default ; https://github.com/CGATOxford/CGATPipelines/pull/422
* update conda from 4.3 to 4.5 (solving "CXXABI_1.3.9' not found" error ; https://github.com/ContinuumIO/anaconda-issues/issues/5191) ; https://github.com/CGATOxford/CGATPipelines/compare/ebc3b7478b774...3ef14607e91b7
* added test for C/C++ compiler ; https://github.com/CGATOxford/CGATPipelines/pull/425

0.3.2
-----

* removed GCProfile calls from pipelines; https://github.com/CGATOxford/CGATPipelines/pull/368
* updated installation; https://github.com/CGATOxford/CGATPipelines/pull/370, https://github.com/CGATOxford/CGATPipelines/pull/371, https://github.com/CGATOxford/CGATPipelines/pull/372, https://github.com/CGATOxford/CGATPipelines/pull/375, https://github.com/CGATOxford/CGATPipelines/pull/377
* pipeline_readqc.py: call cgat script with cgat command; https://github.com/CGATOxford/CGATPipelines/pull/374
* PipelineMapping.py: fix so star will work with muti QC; https://github.com/CGATOxford/CGATPipelines/pull/376
* code changes for improved portability of pipelines; https://github.com/CGATOxford/CGATPipelines/pull/380
* create Python 2 environments for legacy dependencies; https://github.com/CGATOxford/CGATPipelines/pull/381
* re-factored pipeline_annotations and pipeline_mappinginto pipeline_genesets and pipeline_bamstats respectively; https://github.com/CGATOxford/CGATPipelines/pull/373

0.3.1
-----

* fixes for building the docs
* pipeline_peakcalling: changed default blacklist path; https://github.com/CGATOxford/CGATPipelines/pull/345
* pipeline_peakcalling: bedGraphToBigWig in new environment requires sorted input; https://github.com/CGATOxford/CGATPipelines/pull/346
* pipeline_peakcalling: debugging for python 3 in peakcalling notebook reports and changing default format to python 3; https://github.com/CGATOxford/CGATPipelines/pull/347
* added scripts to help find R and Python dependencies: scripts/cgat_deps_R.sh, scripts/cgat_deps_python.sh
* pipeline_peakcalling: added an additional plotting function to generate summary plots after IDR; https://github.com/CGATOxford/CGATPipelines/pull/348
* pipeline_peakcalling: redirected optimal peaks and conservative peaks to filtered IDR output; https://github.com/CGATOxford/CGATPipelines/pull/349
* pipeline_rnaseqqc: added functions to check the strandedness of the library, which is required as an input for various mappers; https://github.com/CGATOxford/CGATPipelines/pull/351
* remove unused imports with autoflakes; https://github.com/CGATOxford/CGATPipelines/pull/354
* bugfix the way bowtie2 index folder is configured in pipeline_mapping; https://github.com/CGATOxford/CGATPipelines/pull/356 and https://github.com/CGATOxford/CGATPipelines/pull/361
* bowtie should not use the quite option so multiqc reports properly; https://github.com/CGATOxford/CGATPipelines/pull/357
* pipeline_peakcalling: minor bugfix; https://github.com/CGATOxford/CGATPipelines/pull/358
* pipeline_peakcalling: update to work with sicer in a Python 2 environment ; https://github.com/CGATOxford/CGATPipelines/pull/359
* added **cgatflow** command; https://github.com/CGATOxford/CGATPipelines/pull/360
* updated documentation and moved unused code to obsolete folder; https://github.com/CGATOxford/CGATPipelines/pull/362
* updated installation; https://github.com/CGATOxford/CGATPipelines/pull/363 ; https://github.com/CGATOxford/CGATPipelines/pull/367
* pipeline_peakcalling: update the way IDR is called so it is included as a conda dependency ; https://github.com/CGATOxford/CGATPipelines/pull/364

0.3.0
-----

* First release working in Python 3 only.

0.1.1
-----

* Last release compatible with Python 2.7

0.1
---

* Re-organize Pipeline.py as a package.
* Move legacy pipelines OPTIC and GPIPE into separate repository.
* Distribute contents of PipelineUtilities to CGAT code collection.
* Revise documentation for modules and scripts.

Contributions
=============

We included publicly and freely available code into the tool
collection for convenience. 

* IGV.py was written by Brent Pedersen.
* SVGdraw.py was written by ...
* The NCL module draws from code written by ...
* list_overlap.py
* Iterators.py

Contributors
============

Andreas Heger
Antonio Berlanga-Taylor
Martin Dienstbier
Nicholas Ilott
Jethro Johnson
Katherine Fawcett
Stephen Sansom
David Sims
Ian Sudbery
Hu Xiaoming
Thomas Smith
Michael Morgan
Katherine Brown
Charlotte George
Adam Cribbs
Hania Pavlou
Reshma Nibhani
Jakub Scaber
Sebastian Luna-Valero

