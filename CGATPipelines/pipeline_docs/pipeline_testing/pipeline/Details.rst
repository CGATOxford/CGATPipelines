.. _filelists:

==========
File lists
==========

The table belows lists for each test the files
that are different, missing or extra.

Different files
===============

Checksums
---------

.. report:: TestingReport.FilesWithProblems
   :render: table
   :groupby: none
   :slices: different_md5

   List of files that have different checksums.

Lines
-----

.. report:: TestingReport.FilesWithProblems
   :render: table
   :groupby: none
   :slices: different_lines

   List of files that have different number of lines


Missing files
=============

.. report:: TestingReport.FilesWithProblems
   :render: table
   :groupby: none
   :slices: missing

   List of missing files in each test.

Extra files
===========

.. report:: TestingReport.FilesWithProblems
   :render: table
   :groupby: none
   :slices: extra

   Extra files in each test.

