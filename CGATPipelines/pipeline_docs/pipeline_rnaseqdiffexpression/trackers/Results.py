import pandas as pd
import numpy as np
import sqlite3

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *

###############################################################################
# parse params
###############################################################################
DATABASE = P.get('', P.get('sql_backend', 'sqlite:///./csvdb'))
ANNOTATIONS_DATABASE = P.get('annotations_database')
###############################################################################
# trackers
###############################################################################


class DeseqFeatureResultsGenes(IsoformTracker):

    pattern = "deseq2_featurecounts__(.*)_genes_results"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.control_name, A.treatment_name, A.control_mean,
        A.treatment_mean, A.test_id, A.l2fold, A.p_value, A.p_value_adj,
        A.significant FROM deseq2_featurecounts__%(track)s_genes_results
        AS A ORDER BY A.significant DESC,
        A.l2fold ASC;
        '''

        return self.getAll(statement)


class EdgerFeatureResultsGenes(IsoformTracker):

    pattern = "edger_featurecounts__(.*)_genes_results"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.control_name, A.treatment_name, A.control_mean,
        A.treatment_mean, A.test_id, A.l2fold, A.p_value, A.p_value_adj,
        A.significant FROM edger_featurecounts__%(track)s_genes_results
        AS A ORDER BY A.significant DESC,
        A.l2fold ASC
        '''

        return self.getAll(statement)


class DeseqKallistoResultsGenes(IsoformTracker):

    pattern = "deseq2_kallisto__(.*)_genes_results"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.control_name, A.treatment_name, A.control_mean,
        A.treatment_mean, A.test_id, A.l2fold, A.p_value, A.p_value_adj,
        A.significant FROM deseq2_kallisto__%(track)s_genes_results
        AS A ORDER BY A.significant DESC,
        A.l2fold ASC
        '''

        return self.getAll(statement)


class EdgerKallistoResultsGenes(IsoformTracker):

    pattern = "edger_kallisto__(.*)_genes_results"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.control_name, A.treatment_name, A.control_mean,
        A.treatment_mean, A.test_id, A.l2fold, A.p_value, A.p_value_adj,
        A.significant FROM edger_kallisto__%(track)s_genes_results AS
        A ORDER BY A.significant DESC,
        A.l2fold ASC
        '''

        return self.getAll(statement)


class SleuthKallistoResultsGenes(IsoformTracker):

    pattern = "sleuth_kallisto__(.*)_genes_results"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.control_name, A.treatment_name, A.control_mean,
        A.treatment_mean, A.test_id, A.l2fold, A.p_value, A.p_value_adj,
        A.significant FROM sleuth_kallisto__%(track)s_genes_results AS
        A ORDER BY A.significant DESC,
        A.l2fold ASC
        '''

        return self.getAll(statement)
