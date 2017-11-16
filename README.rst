.. image:: https://travis-ci.org/CGATOxford/CGATPipelines.svg?branch=master
    :target: https://travis-ci.org/CGATOxford/CGATPipelines

==================
The CGAT pipelines
==================

In CGAT_ we have developed a set of ruffus_ based pipelines in comparative genomics
and NGS analysis. Some documentation of the pipelines is
`here <https://www.cgat.org/downloads/public/cgatpipelines/documentation>`_.

We are working on improving the existing documentation and portability of the code
to release a set of production pipelines soon so please stay tuned.

We are currently testing a script to automate the installation with conda_. Feel
free to give it a go::

        # download installation script:
        curl -O https://raw.githubusercontent.com/CGATOxford/CGATPipelines/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install the development version (recommended, no production version yet):
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # enable the conda environment as requested by the installation script:
        source </full/path/to/folder/without/trailing/slash>/conda-install/bin/activate cgat-p

        # and uninstall pika, which we use internally for our dashboard:
        conda remove pika

        # finally, please run the cgatflow command-line tool to check the installation:
        cgatflow --help

The installation script will put everything under the specified location. It needs
15 GB of disk space and it takes about 35 minutes to complete. The aim of the
script is to provide a portable installation that does not interfere with the existing
software. As a result, you will have a conda environment working with the CGAT Pipelines
which can be enabled on demand according to your needs.

On top of the instructions above, please make sure that you configure the following
environment variables::

        # Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
        export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so

        # You can get this value from your configured environment:
        env | grep DRMAA_LIBRARY_PATH

        # or just look for the library:
        find <path-to-DRMS-install-folder> -name "*libdrmaa.so"

        # Also, make sure you have defined temporary folders
        # 1. Local to execution hosts with
        export TMPDIR=/tmp
        # 2. Shared to pipeline working directory
        export SHARED_TMPDIR=/<path-to-network-folder>/scratch

For questions, please open a new issue on
`GitHub
<https://github.com/CGATOxford/CGATPipelines/issues>`_.

.. _ruffus: http://www.ruffus.org.uk
.. _CGAT: http://www.cgat.org
.. _conda: https://conda.io

