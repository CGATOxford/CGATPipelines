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
class replicatedUniqueIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_replicated_unique_intervals" % locals())
        return odict(list(zip(("Unique intervals", "mean_interval_length"), data)))

##########################################################################


class replicatedUniqueIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT (stop-start) FROM %(track)s_replicated_unique_intervals" % locals())
        return {"length": data}

##########################################################################


class replicatedUniqueIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_replicated_unique_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"peakval": data}

##########################################################################


class replicatedUniqueIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_replicated_unique_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"avgval": data}

##########################################################################


class replicatedUniqueIntervalFoldChange(cpgTracker):

    """Distribution of fold change """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT i.fold FROM %(track)s_replicated_unique_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return odict([("Fold Change", data)])

##########################################################################


class replicatedUniqueIntervalTSS(cpgTracker):

    """Distribution of distance to closest TSS """

    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( '''SELECT closest_dist FROM %(track)s_replicated_unique_intervals u, 
                                  %(track)s_replicated_intervals i, %(track)s_replicated_%(ANNOTATIONS_NAME)s_transcript_tss_distance t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start
                                  AND t.gene_id=i.interval_id''' % locals() )
        return {"distance": data}

##########################################################################


class replicatedUniqueIntervalCpGDensity(cpgTracker):
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT pCpG FROM %(track)s_replicated_unique_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedUniqueIntervalCpGObsExp(cpgTracker):
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT CpG_ObsExp FROM %(track)s_replicated_unique_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedUniqueIntervalCpGNumber(cpgTracker):
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_replicated_unique_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedUniqueIntervalGCContent(cpgTracker):
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_replicated_unique_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedUniqueIntervalTranscriptOverlap(featureOverlap):

    """return overlap of interval with  protein-coding transcripts """
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN  tss_transcript_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap o, %(track)s_replicated_unique_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )

        return odict(list(zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data)))

##########################################################################


class replicatedUniqueIntervalGeneOverlap(featureOverlap):

    """return overlap of interval with  protein-coding genes """
    mPattern = "_replicated_unique_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN tss_gene_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap o, %(track)s_replicated_unique_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )

        return odict(list(zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data)))
