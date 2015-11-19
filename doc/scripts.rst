.. _scripts:

=====================
Scripts
=====================

The script collection of :term:`CGATPipelines` contain scripts that
are useful for setting up and managing pipelines.

Setting up a pipeline
=====================

:doc:`scripts/pipeline_quickstart`
    Create a directory layout for building a new pipeline from
    scratch.

:doc:`scripts/cgat_ruffus_profile`
    Collect benchmarking information from a completed, aborted
    or running pipelines.

:doc:`scripts/cgat_tsv2links`
    Link input files into a pipeline directory for processing and
    relabel them according to the experimental design.

Analysing a pipeline
====================

:doc:`scripts/cgat_logfiles2tsv`
    Collect benchmarking information from job log files.

:doc:`scripts/cgat_build_report_page`
    Build a web-page for all reports build across all CGAT projects.

Cleaning up a pipeline
======================

:doc:`scripts/cgat_clean`
    Remove files from aborted runs in a pipeline directory.

:doc:`scripts/cgat_zap`
    Preserve disk space by zapping intermediate files while at the
    same time preserving pipeline state.

:doc:`scripts/cgat_cwd2list`
    Save a time-stamped list of all files in a pipeline directory.

Cluster and job management
==========================

:doc:`scripts/cgat_cluster_distribute`
    Distribute files on the cluster

:doc:`scripts/run_function`
    Run a function inside a python module on the cluster.

:doc:`scripts/farm`
    Split data into chunks and execute a command on each one
    merging the results.

:doc:`scripts/nofarm`
    Equivalent to :doc:`farm`

:doc:`scripts/qkill`
    Kill queued jobs based on a pattern.

:doc:`scripts/run`
    Execute a command adding meta and benchmark information.

:doc:`scripts/submit`
    Submit a list of qsub scripts to the cluster.
