=============
Release notes
=============

The CGAT pipelines have not been published yet. Please use the latest version from the repository.

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

