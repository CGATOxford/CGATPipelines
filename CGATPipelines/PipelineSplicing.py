
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


def runRMATS(gtffile, designfile, pvalue, strand, outdir, permute=0):
    '''DEExperiment object to generate differential splicing events
    using rMATS
    '''
    design = Expression.ExperimentalDesign(designfile)
    if permute == 1:
        design.table.group = random.choice(list(
                             itertools.permutations(design.table.group)))

    group1 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[0])])
    with open(outdir + "/b1.txt", "w") as f:
        f.write(group1)
    group2 = ",".join(
        ["%s.bam" % x for x in design.getSamplesInGroup(design.groups[1])])
    with open(outdir + "/b2.txt", "w") as f:
        f.write(group1)
    readlength = BamTools.estimateTagSize(design.samples[0]+".bam")

    statement = '''rMATS
    --b1 %(outdir)s/b1.txt
    --b2 %(outdir)s/b2.txt
    --gtf <(gunzip -c %(gtffile)s)
    --od %(outdir)s
    --readLength %(readlength)s
    --cstat %(pvalue)s
    --libType %(strand)s
    ''' % locals()

    # Specify paired design
    if design.has_pairs:
        statement += '''-analysis P '''

    # if Paired End Reads
    if BamTools.isPaired(design.samples[0]+".bam"):
        statement += '''-t paired''' % locals()

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
