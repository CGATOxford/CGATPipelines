'''
======================================================================

Requirements:


Reference
---------

'''

import os
import re
import collections
import pandas
import math
import numpy
import numpy.ma as ma
import itertools
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.IOTools as IOTools

##############################################
# Preprocessing Functions
def filterBams(infile, outfile, filters):

    F = ""
    if 'unmapped' in filters:
        F += str(4)
    if 'unpaired' in filters:
        F += str(12)

    T1 = P.getTempFilename(".")
    T2 = P.getTempFilename(".")
    # ensure bamfile is sorted,
    statement = """samtools sort %(infile)s %(T1)s;
           samtools index %(T1)s;
           samtools view %(T1)s -b -F %(F)s > %(T2)s;"""

    # remove duplicates, if requested
    if 'duplicates' in filters:
        statement += """
        MarkDuplicates \
        INPUT=%(T2)s \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        OUTPUT=%(outfile)s \
        METRICS_FILE=/dev/null \
        VALIDATION_STRINGENCY=SILENT \
        2> %(outfile)s.log"""

    P.run()

#############################################
# Peakcalling Functions

#############################################
# QC Functions
