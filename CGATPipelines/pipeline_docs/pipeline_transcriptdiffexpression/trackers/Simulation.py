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

    def __call__(self, track, slice=None):

        statement = '''
        SELECT read_count, est_counts, fraction_bin, cor
        FROM simulation_correlations
        '''

        return self.getAll(statement)
