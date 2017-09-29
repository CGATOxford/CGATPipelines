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
import cpgReport

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict

##########################################################################
##########################################################################


class intergenicSummary(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_intergenic_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(gene_id) as Intergenic_intervals
                   FROM %(track)s_replicated_intergenic_noncoding'''
        data = self.getAll(query)
        return data

##########################################################################


class noncodingOverlap(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.is_overlap = 1
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data


##########################################################################
class noncoding1kbDist(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.closest_dist < 1000
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data

##########################################################################


class noncoding5kbDist(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.closest_dist < 5000
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data


##########################################################################
class noncodingOverlapIntergenic(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_intergenic_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_intergenic_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.is_overlap = 1
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data


##########################################################################
class noncoding1kbDistIntergenic(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_intergenic_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_intergenic_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.closest_dist < 1000
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data

##########################################################################


class noncoding5kbDistIntergenic(cpgReport.cpgTracker):

    '''Summary table'''
    mPattern = "_replicated_intergenic_noncoding$"

    def __call__(self, track, slice=None):
        query = '''SELECT count(n.gene_id) as Intervals, t.gene_biotype
                   FROM %(track)s_replicated_intergenic_noncoding n, annotations.transcript_info t
                   WHERE substr(n.closest_id,1,18)=t.gene_id
                   AND n.closest_dist < 5000
                   GROUP BY gene_biotype
                   ORDER BY Intervals desc'''
        data = self.getAll(query)
        return data


##########################################################################
##########################################################################
##########################################################################
class noncodingTSSOverlap(cpgReport.cpgTracker):

    '''number of transcript TSSs that an interval overlaps.'''
    mPattern = "_replicated_noncoding$"
    mAnnotations = "replicated_annotations"
    mTable = "replicated_noncoding"
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


class noncodingTSSClosest(cpgReport.cpgTracker):

    """for each interval, return the distance to the closest TSS."""

    mXLabel = "distance / bases"
    mPattern = "_replicated_noncoding$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "replicated_annotations"
    mTable = "replicated_noncoding"

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
