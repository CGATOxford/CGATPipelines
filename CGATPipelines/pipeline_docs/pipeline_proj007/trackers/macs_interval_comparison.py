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


class ExternalIntervalLists(cpgTracker):

    """Summary stats of external interval lists used for comparison. """

    mPattern = "external_interval_sets$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT bed, intervals FROM external_interval_sets" % locals())
        return odict(list(zip(("Dataset", "Intervals"), list(zip(*data)))))

##########################################################################


class OverlapCpG(cpgTracker):

    """Count of intervals overlapping CpG islands for each dataset. """

    mPattern = "_replicated_cgi_venn$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_replicated_cgi_venn c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_replicated_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(list(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), list(zip(*data)))))

##########################################################################


class OverlapChipseq(cpgTracker):

    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_replicated_chipseq$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track,b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_replicated_chipseq c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_replicated_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(list(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), list(zip(*data)))))

##########################################################################


class OverlapCAPseq(cpgTracker):

    """Count of intervals overlapping CAPseq intervals for each dataset. """

    mPattern = "_replicated_capseq$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_replicated_capseq c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_replicated_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(list(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), list(zip(*data)))))

##########################################################################


class OverlapChromatinMarks(cpgTracker):

    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_replicated_chromatin$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_replicated_chromatin c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_replicated_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(list(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), list(zip(*data)))))

##########################################################################


class gatResults(cpgTracker):

    """Summary stats of GAT analysis. """

    mPattern = "gat_results$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM external_dataset_gat_results ")
        return odict(list(zip(("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value"), list(zip(*data)))))
