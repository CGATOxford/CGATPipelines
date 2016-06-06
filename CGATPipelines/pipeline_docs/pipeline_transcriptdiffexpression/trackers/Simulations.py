import pandas as pd
import numpy as np

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *


class imagesTracker(TrackerImages):

    '''Convience Tracker for globbing images for gallery plot'''
    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]


class simulationCorrelations(IsoformTracker):

    pattern = "(\S+)_simulation_correlations$"
    metric = ""

    def __call__(self, track, slice=None):

        statement = '''
        SELECT tpm, read_count, est_%(metric)s, fraction_unique,
        fraction_bin, %(metric)s_cor,
        log2diff_%(metric)s, log2diff_%(metric)s_thres
        FROM %(track)s_simulation_correlations
        ORDER BY log2diff_%(metric)s
        '''

        return self.getAll(statement)


class simulationCorrelationsCount(simulationCorrelations):

    pattern = "(\S+)_simulation_correlations$"
    metric = "counts"


class simulationCorrelationsTpm(simulationCorrelations):

    pattern = "(\S+)_simulation_correlations$"
    metric = "tpm"


class simulationCorrelationsSummaryFold(IsoformTracker):

    pattern = "(\S+)_simulation_correlations$"
    metric = ""

    def __call__(self, track, slice=None):

        statement = '''
        SELECT
        total(CASE WHEN abs(log2diff_%(metric)s) >1 THEN 1 ELSE 0 END)
        AS 'flagged',
        total(CASE WHEN abs(log2diff_%(metric)s) <1 THEN 1 ELSE 0 END)
        AS 'passed'
        FROM %(track)s_simulation_correlations
        '''

        return self.getAll(statement)


class simulationCorrelationsSummaryFoldCount(simulationCorrelationsSummaryFold):
    pattern = "(\S+)_simulation_correlations$"
    metric = "counts"


class simulationCorrelationsSummaryFoldTpm(simulationCorrelationsSummaryFold):
    pattern = "(\S+)_simulation_correlations$"
    metric = "tpm"


class simulationCorrelationsSummaryKmers(simulationCorrelationsSummaryFold):

    pattern = "(\S+)_simulation_correlations$"

    def __call__(self, track, slice=None):

        statement = '''
        select
        total(CASE WHEN fraction_unique <0.03 THEN 1 ELSE 0 END) AS 'flaged',
        total(CASE WHEN fraction_unique >=0.03 THEN 1 ELSE 0 END) AS 'passed'
        FROM %(track)s_simulation_correlations
        '''

        return self.getAll(statement)
