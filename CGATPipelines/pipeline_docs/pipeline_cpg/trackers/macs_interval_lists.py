import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import scipy.stats
import numpy.ma
import Stats
import Histogram

from CGATReport.Tracker import *
from cpgReport import *

##########################################################################


class IntervalList(cpgTracker):

    '''list of intervals.'''

    nresults = 20
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold")
    mPattern = "-?_macs_intervals$"

    def getSQLStatement(self, track, slice=None):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY i.avgval DESC''' % locals()

        if self.nresults:
            statement += " LIMIT %i" % self.nresults

        return statement

    def __call__(self, track, slice=None):

        statement = self.getSQLStatement(track, slice)
        data = self.get(statement)
        ucsc_genome = UCSC_GENOME
        n = odict()
        for d in data:
            id, contig, start, end, length = d[:5]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_genome)s&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[str(id)] = odict(
                list(zip(self.mColumnsFixed + self.mColumnsVariable, (pos, length,) + d[5:])))

        return n

##########################################################################


class IntervalListFull(cpgTracker):

    '''list of all intervals. Table for export. '''

    nresults = None
    mPattern = "_macs_intervals$"

    def __call__(self, track, slice=None):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY i.peakval DESC''' % locals()

        data = self.get(statement)
        return odict(list(zip(("contig", "start", "end", "peakval", "avgval"),  list(zip(*data)))))

##########################################################################


class IntervalListPeakval(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY i.peakval DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListAvgval(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY i.avgval DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListFoldChange(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY fold DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListLength(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_intervals AS i
                       ORDER BY length DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class LongIntervals(IntervalList):

    '''list of intervals >10kb in length sorted by fold change'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, round(i.avgval,2), i.fold
                       FROM %(track)s_macs_intervals AS i
                       WHERE i.length > 10000
                       ORDER BY fold DESC''' % locals()
        return statement

##########################################################################


class IntervalListLowGC(IntervalList):

    '''list of intervals with pGC <0.5 and pCpG <0.6 sorted by fold change'''

    nresults = 100
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold", "pGC", "CpG_ObsExp2")

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, round(i.avgval,2), i.fold, 
                       round(c.pGC,2), round(c.CpG_ObsExp2,2)
                       FROM %(track)s_macs_intervals AS i, %(track)s_composition AS c
                       WHERE i.interval_id=c.gene_id
                       AND c.CpG_ObsExp2 < 0.6
                       AND c.pGC < 0.5
                       ORDER BY i.fold DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListCDS(IntervalList):

    '''list of intervals overlapping CDS.'''
    nresults = 100
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold", "nover_CDS",
                        "pover1_CDS", "pover2_CDS", "closest_id", "gene_id", "gene_name")
    mPattern = "-?_macs_intervals$"

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT distinct i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold, 
                       a.nover_CDS, a.pover1_CDS, a.pover2_CDS, 
                       tss.closest_id, tr.gene_id, tr.gene_name
                       FROM %(track)s_macs_intervals AS i, %(track)s_annotations AS a,
                       %(track)s_tss AS tss, annotations.transcript_info AS tr
                       WHERE i.interval_id=a.gene_id 
                       AND i.interval_id=tss.gene_id
                       AND tr.transcript_id=tss.closest_id
                       AND a.nover_CDS>0
                       ORDER BY a.pover2_CDS DESC, a.pover1_CDS DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class FoldChangeCounts(cpgTracker):

    """Correlation between all sets. """

    pattern = "(.*)_macs_intervals$"
    mMinFoldChange = 2.0

    def __call__(self, track, slice=None):
        data = []
        upfold = self.mMinFoldChange
        downfold = 1.0 / upfold
        data.append(("> %5.2f fold" % upfold, self.getValue(
            "SELECT COUNT(*) FROM %(track)s_macs_intervals WHERE fold > %(upfold)f " % locals())))
        data.append(("unchanged", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_macs_intervals WHERE fold between %(downfold)f and %(upfold)f" % locals())))
        data.append(("< %5.2f fold" % downfold, self.getValue(
            "SELECT COUNT(*) FROM %(track)s_macs_intervals WHERE fold < %(downfold)f " % locals())))

        return odict(data)
