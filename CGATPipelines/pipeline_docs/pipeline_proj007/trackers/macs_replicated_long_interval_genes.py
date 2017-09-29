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
from macs_interval_lists import *
from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict

##########################################################################
##########################################################################
##########################################################################
##########################################################################


class LongIntervals(IntervalList):

    '''list of intervals >3kb in length sorted by fold change'''

    nresults = 1000
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = (
        "avgval", "fold", "gene_id", "gene_name", "genes_pover1", "genes_pover2")
    mPattern = "_replicated_intervals$"

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults
        ANNOTATIONS_NAME = P['annotations_name']
        statement = '''SELECT distinct i.interval_id, i.contig, i.start, i.end, i.length, i.avgval, i.fold, t.gene_id, t.gene_name, o.genes_pover1, o.genes_pover2
                       FROM %(track)s_replicated_intervals i, %(track)s_replicated_%(ANNOTATIONS_NAME)s_transcript_tss_distance s, 
                       %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_transcript_mapping m,
                       %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap o, annotations.transcript_info t
                       WHERE i.interval_id=s.gene_id
                       AND o.gene_id=i.interval_id
                       AND i.length > 3000
                       AND o.genes_pover2 > 0
                       AND t.gene_biotype='protein_coding'
                       AND s.gene_id=m.interval_id
                       AND m.transcript_id=t.transcript_id
                       ORDER BY o.length desc
                       LIMIT 500''' % locals()
        print(statement)
        return statement

##########################################################################


class longGenesCapseqProfile(TrackerImages):

    """CAPseq profile per gene """

##########################################################################


class longGenesH3K27Profile(TrackerImages):

    """Chromatin profile per gene"""

##########################################################################


class longGenesH3K4Profile(TrackerImages):

    """Chromatin profile per gene"""

##########################################################################


class shortGenesH3K27Venn(TrackerImages):

    """intersection of short genes and H3K27Me3 intervals """

##########################################################################


class longGenesH3K27Venn(TrackerImages):

    """intersection of long genes and H3K27Me3 intervals"""

##########################################################################


class longGenesTissueVenn(TrackerImages):

    """Conservation of long genes across tissues"""

##########################################################################


class longPolycombGAT(cpgTracker):

    """genomic assocation of H3K27Me3 intervals and genes long >90% by NMIs"""
    mPattern = "long_intervals_gat_results$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM long_intervals_gat_results ")
        return odict(list(zip(("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value"), list(zip(*data)))))

##########################################################################


class longPolycombIntersection(cpgTracker):

    """Intersection of H3K27Me3 intervals and genes overlapped by NMIs >3kb in length"""
    mPattern = "long_interval_h3k27me3_venn$"

    def __call__(self, track, slice=None):
        query = '''select * from long_interval_h3k27me3_venn'''
        data = self.getAll(query)
        return data

##########################################################################


class longGenesGOAnalysisBP(cpgTracker):

    '''GO analysis biological process'''
    mPattern = "_long_go_biol_process$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, spercent as percent_of_list, fdr
                   from %(track)s_long_go_biol_process
                   where fdr < 0.05
                   order by fdr asc, scount desc'''
        data = self.getAll(query)
        return data

##########################################################################


class longGenesGOAnalysisCL(cpgTracker):

    '''GO analysis cell location'''
    mPattern = "_long_go_cell_location$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, spercent as percent_of_list, fdr 
                   from %(track)s_long_go_cell_location
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class longGenesGOAnalysisMF(cpgTracker):

    '''GO analysis molecular function'''
    mPattern = "_long_go_mol_function$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, spercent as percent_of_list, fdr 
                   from %(track)s_long_go_mol_function
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class longGenesGOSlimAnalysisBP(cpgTracker):

    '''GO slim analysis biological process '''
    mPattern = "_long_goslim_biol_process$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, spercent as percent_of_list, fdr
                   from %(track)s_long_goslim_biol_process
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data
