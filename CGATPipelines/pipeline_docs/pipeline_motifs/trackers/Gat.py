import os
import sys
import re
import types
import itertools
import math

from IntervalReport import *


class GatTracker(IntervalTracker):
    pass

# class GatGenomicContextTable( GatTracker ):
#     pattern = "gat_context_(.*)$"

#     def __call__(self, track):
#         return self.getAll( "SELECT * FROM gat_context_%(track)s" )

# class GatGenomicAnnotationTable( GatTracker ):
#     pattern = "gat_annotations_(.*)$"

#     def __call__(self, track):
#         return self.getAll( "SELECT * FROM gat_annotations_%(track)s" )

##########################################################################
##########################################################################
##########################################################################
# GAT results
##########################################################################
# class GatResults( IntervalTracker, SingleTableTrackerRows ):
#     '''All gat results.'''
#     fields = ('track', 'annotation')
#     extra_columns = { "colour" : "CASE WHEN qvalue < 0.05 THEN 'red' ELSE 'blue' END" }
#     sort = 'l2fold'


class GatFold(IntervalTracker, SingleTableTrackerEdgeList):

    '''fold change matrix.'''
    row = "track"
    column = "annotation"
    value = "fold"
    where = "pvalue < 0.05"


class GatLogFold(IntervalTracker):

    '''logfold - colour is signficance'''
    fdr = 2.0
    as_tables = True

    def __call__(self, track):
        return self.getDict( """SELECT annotation, fold, 
                                       CASE WHEN qvalue < %(fdr)f THEN 'red' ELSE 'blue' END AS colour
                               FROM %(track)s ORDER BY fold""")


class GatResults(GatTracker):
    as_tables = True

    def __call__(self, track):
        return self.getAll("SELECT * FROM %(track)s")


class GatTableAnnotations:
    pattern = "gat_annotations_(.*)"


class GatTableContext:
    pattern = "gat_context_(.*)"


class GatTableFunctions:
    pattern = "gat_functions_(.*)"

_gat_analysis = {"Results": GatResults,
                 "Fold": GatLogFold,
                 "LogFold": GatLogFold}

_gat_sets = {"Annotations": GatTableAnnotations,
             "Context": GatTableContext,
             "Functions": GatTableFunctions,
             }

for a, aa in list(_gat_analysis.items()):
    for b, bb in list(_gat_sets.items()):
        n = "Gat%s%s" % (a, b)
        globals()[n] = type(n, (bb, aa), {})
