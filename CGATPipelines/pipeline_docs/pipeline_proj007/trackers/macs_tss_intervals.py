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


class tssIntervalSummary(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getFirstRow("SELECT COUNT(t.gene_id) as tss FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                                 ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals())
        data2 = self.getFirstRow("SELECT COUNT(t.gene_id) as tss FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                                 ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals())
        return odict(list(zip(("TSS intervals", "Non TSS intervals"), data1 + data2)))

##########################################################################


class tssIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT i.end-i.start as length FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT i.end-i.start as length FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}

##########################################################################


class tssIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT i.peakval FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT i.peakval FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}

##########################################################################


class tssIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT i.avgval FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT i.avgval FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}

##########################################################################


class tssIntervalFoldChange(cpgTracker):

    """Distribution of fold change """

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT i.fold FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT i.fold FROM %(track)s_replicated_intervals i, %(track)s_replicated_" +
                               ANNOTATIONS_NAME + "_transcript_tss_distance t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals())
        # return odict( [("TSS interval length", data1),("Non-TSS interval
        # length",data2)] )
        return {"TSS interval": data1, "Non-TSS interval": data2}


##########################################################################
class tssIntervalCpGDensity(cpgTracker):

    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT c.pCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT c.pCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}

##########################################################################


class tssIntervalCpGObsExp2(cpgTracker):
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT c.CpG_ObsExp FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT c.CpG_ObsExp FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals())
        return {"TSS interval": data1, "Non-TSS interval": data2}

##########################################################################


class tssIntervalCpGNumber(cpgTracker):
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT c.nCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT c.nCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distances t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}

##########################################################################


class tssIntervalGCContent(cpgTracker):
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        data1 = self.getValues("SELECT c.pGC FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals())
        data2 = self.getValues("SELECT c.pGC FROM %(track)s_replicated_intervals i, %(track)s_replicated_" + ANNOTATIONS_NAME +
                               "_transcript_tss_distance t, %(track)s_replicated_capseq_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals())
        return {"TSS-associated CAPseq intervals": data1, "Non-TSS associated CAPseq intervals": data2}
