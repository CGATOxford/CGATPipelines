.. _status:

======
Status
======

MD5 Comparison
==============

The MD5 Comparison compares the files output by the latest pipeline
run against the results of a reference pipeline run. A run was
successfull if exactly the same files are present and all files are
identical. A list of the files are found in the :ref:`filelists`.

.. report:: Status.ComparisonStatus
   :render: status

   Pipeline Status

Completion status
====================

Looking at the logfiles, the this section lists the completion status
of runnning the pipeline and building the report.

.. report:: Status.PipelineStatus 
   :render: status                  
   
   Status report of pipeline running. The info field shows the last
   time the pipeline was started.

Report status
=============

Looking at the logfiles of the report building, the report below
list the number of errors and warnings encountered when building the
report.

.. report:: TestingReport.ReportTable
   :render: table

   Table summarizing the logfile of running the reports. The columns
   are:
   report_error -- number of errors in the report
   report_warning -- number of warnings in the report
   error -- number of error messsages in the logfile
   warning -- number of warning messages in the logfile
