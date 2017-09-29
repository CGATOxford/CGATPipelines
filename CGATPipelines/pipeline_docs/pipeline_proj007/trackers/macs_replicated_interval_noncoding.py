import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
from cpgReport import *
from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict

##########################################################################
##########################################################################


class noncodingOverlap(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.is_overlap = 1
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        # print query
        data = self.getAll(query)
        return data

##########################################################################


class noncoding1kbDist(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.closest_dist < 1000
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data

##########################################################################


class noncoding5kbDist(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.closest_dist < 5000
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data

##########################################################################
##########################################################################
##########################################################################


class lncrnaCount(cpgTracker):

    '''Summary table'''
    mPattern = "lncrna_bed$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct id) as lncRNAs FROM lncrna_bed'''
        data = self.getAll(query)
        return data

##########################################################################


class lncrnaOverlap(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_lncrna_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, count(distinct m.gene_id) as lncRNAs
                   FROM %(track)s_replicated_lncrna_tss_distance n, %(track)s_replicated_interval_lncrna_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.is_overlap = 1'''
        print(query)
        data = self.getAll(query)
        return data

##########################################################################


class lncrna1kbDist(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_lncrna_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, count(distinct m.gene_id) as lncRNAs
                   FROM %(track)s_replicated_lncrna_tss_distance n, %(track)s_replicated_interval_lncrna_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.closest_dist < 1000'''
        data = self.getAll(query)
        return data

##########################################################################


class lncrna5kbDist(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_lncrna_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as Intervals, count(distinct m.gene_id) as lncRNAs
                   FROM %(track)s_replicated_lncrna_tss_distance n, %(track)s_replicated_interval_lncrna_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.closest_dist < 5000'''
        data = self.getAll(query)
        return data

##########################################################################
##########################################################################
##########################################################################


class noncodingTSSOverlap(cpgTracker):

    '''number of transcript TSSs that an interval overlaps.'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"
    mAnnotations = "_annotations"
    mTable = "replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance"
    mColumn = "d.is_overlap"
    mWhere = "d.is_overlap < 5 "

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        hist, bins = numpy.histogram(
            data, bins=numpy.arange(0, max(data) + 1, 1))
        return odict(list(zip(list(map(str, bins[:-1])), hist)))

##########################################################################


class noncodingTSSClosest(cpgTracker):

    """for each interval, return the distance to the closest TSS."""
    ANNOTATIONS_NAME = P['annotations_name']
    mXLabel = "distance / bases"
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "_annotations"
    mTable = "replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.get(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.get( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        return data

##########################################################################


class noncodingTSSClosestUpstream(noncodingTSSClosest):

    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist5 is null THEN 1000000 ELSE dist5 END as dist5 "

##########################################################################


class noncodingTSSClosestDownstream(noncodingTSSClosest):

    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist3 is null THEN 1000000 ELSE dist3 END as dist3 "

##########################################################################


class noncodingCapseqOverlap(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct t.gene_id) as genes, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.is_overlap = 1
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY genes desc'''
        # print query
        data = self.getAll(query)
        return data

##########################################################################


class noncodingCapseq1kbDist(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct t.gene_id) as genes, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.closest_dist < 1000
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY genes desc'''
        # print query
        data = self.getAll(query)
        return data

##########################################################################


class noncodingCapseq5kbDist(cpgTracker):

    '''Summary table'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_noncoding_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct t.gene_id) as genes, t.gene_biotype
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n, 
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m, 
                   annotations.transcript_info t
                   WHERE m.gene_id=t.gene_id
                   AND m.interval_id=n.gene_id
                   AND n.closest_dist < 5000
                   AND t.gene_biotype != "protein_coding"
                   GROUP BY gene_biotype
                   ORDER BY genes desc'''
        # print query
        data = self.getAll(query)
        return data

##########################################################################
##########################################################################


class rnaseqCount(cpgTracker):

    '''Summary table'''
    mPattern = "rnaseq_bed$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct id) as RNAseq_transcripts FROM rnaseq_bed'''
        data = self.getAll(query)
        return data

##########################################################################


class rnaseqOverlap(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_rnaseq_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as NMIs, count(distinct m.gene_id) as Transcripts
                   FROM %(track)s_replicated_rnaseq_tss_distance n, 
                   %(track)s_replicated_interval_rnaseq_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.is_overlap = 1'''
        # print query
        data = self.getAll(query)
        return data

##########################################################################


class rnaseq1kbDist(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_rnaseq_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as NMIs, count(distinct m.gene_id) as Transcripts
                   FROM %(track)s_replicated_rnaseq_tss_distance n, 
                   %(track)s_replicated_interval_rnaseq_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.closest_dist < 1000'''
        data = self.getAll(query)
        return data

##########################################################################


class rnaseq5kbDist(cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_rnaseq_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(distinct n.gene_id) as NMIs, count(distinct m.gene_id) as Transcripts
                   FROM %(track)s_replicated_rnaseq_tss_distance n, 
                   %(track)s_replicated_interval_rnaseq_mapping m
                   WHERE m.interval_id=n.gene_id
                   AND n.closest_dist < 5000'''
        data = self.getAll(query)
        return data
