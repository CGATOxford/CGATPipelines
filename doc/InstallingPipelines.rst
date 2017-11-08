.. _CGATSetup:

=========================
Installation instructions
=========================

The section below describes how to install the CGAT_ Pipelines. We distinguish between two
different installation types: production and development. The former refers to a well
tested subset of pipelines, and is the recommended installation. The latter refers to
the whole collection of pipelines developed at CGAT, which may contain code under active
development.

Please note that we can not test our code on all systems and configurations out there so
please bear with us.

Automated installation
======================

The preferred method to install the CGAT Pipelines is using the installation script,
which uses conda_.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/CGATOxford/CGATPipelines/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install set of production scripts (well tested):
        bash install-CGAT-tools.sh --production [--location </full/path/to/folder/without/trailing/slash>]

        # or go for the latest development version:
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

The installation script will put everything under the specified location. The aim of the
script is to provide a portable installation that does not interfere with the existing
software. As a result, you will have a conda environment working with the CGAT Pipelines
which can be enabled on demand according to your needs.

Manual installation
===================

To obtain the latest code, check it out from the public git_ repository and activate it::

   git clone https://github.com/CGATOxford/CGATPipelines.git
   cd CGATPipelines
   python setup.py develop

The CGAT Pipelines depends on the CGAT scripts, which can be installed by following the
installation instructions `here
<http://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html>`_.

Once checked-out, you can get the latest changes via pulling::

   git pull 

Some scripts contain cython_ code that needs to be recompiled if the
script or the pysam_ installation has changed. To rebuild all scripts,
for example after updating the repository, type::

   python cgat/scripts/cgat_rebuild_extensions.py

Recompilation requires a C compiler to be installed.

Setting up the computing environment
====================================

The pipelines assume that Sun Grid Engine has been installed. Other
queueing systems might work, but expect to be disappointed. The
pipeline is started on a :term:`submit host` assuming a default queue
``all.q``. Other queues can be specified on the command line, for
example::

    python cgat/CGATPipelines/pipeline_<name>.py --cluster-queue=medium_jobs.q

A pipeline might start up to ``-p/--multiprocess`` processes. Preferentially,
tasks are sent to the cluster, but for some tasks this is not possible. 
These might thus run on the :term:`submit host`, so make sure it is fairly powerful.

Pipelines expects that the :term:`working directory` is accessible with
the same path both from the submit and the :term:`execution host`. 

Software requirements
=====================

CGAT Pipelines make use of a variety of software. We keep a list of software dependencies
in the form of a conda_ environment file `here
<https://github.com/CGATOxford/CGATPipelines/blob/master/conda/environments/pipelines-devel.yml>`_.

All these dependencies will be automatically installed with the automated installation
script as explained above.

What exactly is required will depend on the particular pipeline. The
pipeline assumes that the executables are in the users :envvar:`PATH`
and that the rest of the environment has been set up for each tool.

To check if the dependencies within a particular pipeline are satisfied, type::

   python CGATPipelines/pipeline_mapping.py --input-validation

.. _conda: https://conda.io
.. _CGAT: http://www.cgat.org
