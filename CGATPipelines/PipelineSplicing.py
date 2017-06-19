
#########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
PipelineSplicing.py - wrap various differential expression tools
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides tools for differential splicing analysis.

Methods implemented are:

   rMATS

DEXSeq is implemented elsewhere (counts2table), as it is closely related to
counts-based differential expression tools.

Usage
-----

ru
Documentation
-------------

Requirements:

* rMATS >= ?



'''

from ruffus import *
import os
import random
import itertools
import CGAT.BamTools as BamTools
import CGAT.Experiment as E
import CGAT.Expression as Expression
import CGATPipelines.Pipeline as P


def runRMATS(gtffile, designfile, pvalue, strand, outfile, permute=0):
    '''DEExperiment object to generate differential splicing events
    using rMATS
    '''
    design = Expression.ExperimentalDesign(designfile)
    if permute == 1:
        design.table.group = random.choice(list(
                             itertools.permutations(design.table.group)))

    group1 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[0])])
    group2 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[1])])
    readlength = BamTools.estimateTagSize(design.samples[0]+".bam")
    # outfile = os.path.dirname(os.path.dirname(outfiles[0]))

    statement = '''rMATS
    -b1 %(group1)s
    -b2 %(group2)s
    -gtf <(gunzip -c %(gtffile)s)
    -o %(outfile)s
    -len %(readlength)s
    -c %(pvalue)s
    -libType %(strand)s
    ''' % locals()

    # Specify paired design
    if design.has_pairs:
        statement += '''-analysis P '''

    # Get Insert Size Statistics if Paired End Reads
    if BamTools.isPaired(design.samples[0]+".bam"):
        inserts1 = [BamTools.estimateInsertSizeDistribution(sample+".bam",
                                                            10000)
                    for sample in design.getSamplesInGroup(design.groups[0])]
        inserts2 = [BamTools.estimateInsertSizeDistribution(sample+".bam",
                                                            10000)
                    for sample in design.getSamplesInGroup(design.groups[1])]
        r1 = ",".join(map(str, [item[0] for item in inserts1]))
        sd1 = ",".join(map(str, [item[1] for item in inserts1]))
        r2 = ",".join(map(str, [item[0] for item in inserts2]))
        sd2 = ",".join(map(str, [item[1] for item in inserts2]))

        statement += '''-t paired
        -r1 %(r1)s -r2 %(r2)s
        -sd1 %(sd1)s -sd2 %(sd2)s''' % locals()

    P.run()


def rmats2sashimi(infile, designfile, FDR, outfile):

    Design = Expression.ExperimentalDesign(designfile)
    if len(Design.groups) != 2:
        raise ValueError("Please specify exactly 2 groups per experiment.")

    g1 = Design.getSamplesInGroup(Design.groups[0])
    g2 = Design.getSamplesInGroup(Design.groups[1])

    if len(g1) != len(g2):
        g1 = g1[:min(len(g1), len(g2))]
        g2 = g2[:min(len(g1), len(g2))]
        E.info("The two groups compared were of unequal size. For  " +
               "visual display using sashimi they have been truncated " +
               "to the same length")

    group1 = ",".join(["%s.bam" % x for x in g1])
    group2 = ",".join(["%s.bam" % x for x in g2])
    group1name = Design.groups[0]
    group2name = Design.groups[1]
    event = os.path.basename(os.path.normpath(outfile))

    statement = '''cat
    %(infile)s|grep -v NA|
    awk '$20 < %(FDR)s' > %(infile)s_sig.txt;
    checkpoint;
    rmats2sashimiplot
    --b1 %(group1)s
    --b2 %(group2)s
    -t %(event)s
    -e %(infile)s_sig.txt
    --l1 %(group1name)s
    --l2 %(group2name)s
    -o %(outfile)s; checkpoint;
    ''' % locals()

    P.run()
