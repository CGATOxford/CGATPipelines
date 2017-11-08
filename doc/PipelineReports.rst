.. _WritingReports:

========================
Writing pipeline reports
========================

In CGAT we have a number of ways of generating and implimenting reports in our pipelines.Currently,
we have 4 ways of generating reports:

   * MultiQC
   * Rmarkdown
   * Jupyter notebook
   * CGAT report

Below we have information on how to build and render a report within a pipeline.

Our basic method of rendering a report is to run the `build_report` workflow task::

   cgatflow make build_report

This will impliment an associated report.For our upstream generic pipleines this
is usually multiQC.

**MultiQC**
===========

Writing a report
----------------

In its present format multiQC is not easily adaptable to running bespoke report analysis
tools. Therefore we usually use the default multiQC which navigated through the dictrectory structure and identifies
log files that it recognises. More information about supported tools can be found here: http://multiqc.info/.


Generating a report in the pipeline
-----------------------------------



**Rmarkdown**
=============

Writing a report
----------------

An example of a Rmarkdown report implimentation can be found in
the bamstats pipeline. Detailed instructions on how to genrate a
report can be found on the RStudio website. 


Generating a report in the pipeline
-----------------------------------


**Jupyter Notebook**
====================

Writing a report
----------------

An example of a Jupyter notebook report implimentation can be found in
the bamstats pipeline. 

Generating a report in the pipeline
-----------------------------------



**CGAT Report**
===============

CGAT report is slowly being depricated from a number of our pipelines.
We are hoping to replace all of our reports with Jupyter notebook or
Rmarkdown implimentations.

CGAT pipelines use CGATReport_ to report the outcome of a pipeline
run. Conceptually, the workflow is that a CGAT pipeline creates data
and uploads it into a database. CGATReport_ then creates a report
from the database.


Conditional content
-------------------

The ifconfig_ extension allows to include content depending on configuration
values. To use this extension you will need to modify
:file:`conf.py`. The example below shows the modifications implemented
in :doc:`pipelines/pipeline_mapping` to permit the conditional
inclusion of sections of the report depending on the mapper chosen::

    # add sphinx.ext.ifconfig to the list of extensions
    extensions.append( 'sphinx.ext.ifconfig' )
    
    # define a new configuration variable
    ################################################################
    # Add custom configuration variables for ifconfig extension
    def setup(app):
    	app.add_config_value('MAPPERS', '', True)

    # Set the value of custom configuration variables
    import CGATPipelines.Pipeline as P
    P.getParameters(
	["%s/pipeline.ini" % os.path.splitext(__file__)[0],
	     "../pipeline.ini",
	     "pipeline.ini" ] )

    MAPPERS = P.asList( P.PARAMS["mappers" ] )
    
The thus defined and set custom configuration value ``MAPPERS`` can
now be used inside an rst document::

   .. toctree::
      :maxdepth: 2

      pipeline/Methods.rst
      pipeline/Status.rst
      pipeline/Mapping.rst
      pipeline/MappingSummary.rst
      pipeline/MappingContext.rst
      pipeline/MappingAlignmentStatistics.rst
      pipeline/MappingComplexity.rst

   .. ifconfig:: "tophat" in MAPPERS

      .. toctree::
	 pipeline/MappingTophat.rst

   .. ifconfig:: "star" in MAPPERS

      .. toctree::
	 pipeline/MappingStar.rst

   .. ifconfig:: "tophat" in MAPPERS or "star" in MAPPERS or "gsnap" in MAPPERS

      .. toctree::
	 pipeline/Validation.rst

Note that ``.. ifconfig`` needs to be a first level directive and
can not be include into another directive such as ``.. toctree``.

Referring to other reports
--------------------------

.. note::

   The following below still works, but is obsolete. Instead
   edit your :file:`pipeline.ini` file to include the intersphinx
   mapping targets.


The intersphinx_ extension permits referring to other
CGATReport_ documents. To use this extension you will need to include
the intersphinx_ extension in your :file:`conf.py` configuration file::

    # add sphinx.ext.ifconfig to the list of extensions
    extensions.append( 'sphinx.ext.intersphinx' )

Next, you can add a section called ``intersphinx`` to
:file:`pipeline.ini`::

   [intersphinx]
   readqc=/ifs/projects/proj013/readqc/report/html
   mapping1=/ifs/projects/proj013/mapping1/report/html
   mapping2=/ifs/projects/proj013/mapping2/report/html

.. note::

   It is also possible to add an intersphinx mapping to :file:`conf.py`::

     # add mapping information
     intersphinx_mapping = {
	'readqc': ('/ifs/projects/proj013/readqc/report/html', None) ,
	'mapping1': ('/ifs/projects/proj013/mapping1/report/html', None),
	'mapping2': ('/ifs/projects/proj013/mapping2/report/html', None),
	 }

   The benefit of using :file:`pipeline.ini` is that when a report is
   published :doc:`pipelinemodules.Pipeline` is aware of the links and will
   update the file URLs to web URLs.
   	
This will link to three other reports. The three reports are
abbreviated as ``readqc``, ``mapping1`` and ``mapping2``. The paths
need to be the absolute location of the html build of the sphinx
documents you created previously. These directories should contain a
:file:`objects.inv` file which is usually automatically created by sphinx.

To refer to the other documentation, type::

    :ref:`My link to another documentation <identifier:label>`

``label`` is a valid identifier in the referred to
document. For example::

    :ref:`ReadQC <readqc:readqcpipeline>`

	ReadQC pipeline - fastqc

    :ref:`Unique Mapping  <mapping1:mappingpipeline>`

	Mapping pipeline - short read mapping with bwa. Only
	uniquely mapping reads are kept.

    :ref:`Non-unique mapping <mapping2:mappingpipeline>`

	Mapping pipeline - short read mapping with bwa with same
	parameters as above, but all reads are kept. 

.. _intersphinx: http://sphinx-doc.org/ext/intersphinx.html
.. _ifconfig: http://sphinx-doc.org/ext/ifconfig.html



