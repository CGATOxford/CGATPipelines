import pandas as pd
import numpy as np

from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

from IsoformReport import *

###############################################################################
# parse params
###############################################################################
DATABASE = P.get('', P.get('sql_backend', 'sqlite:///./csvdb'))

###############################################################################
# trackers
###############################################################################

class SleuthResults(IsoformTracker):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.gene_name, A.gene_id, A.transcript_id, B.reason AS flagged,
        A.control_mean AS expression, A.fold, A.l2fold, A.p_value,
        A.p_value_adj, A.significant, A.transcript_biotype
        FROM %(track)s_DEresults AS A
        LEFT JOIN kallisto_flagged_transcripts AS B
        ON A.test_id = B.transcript_id
        %(where)s
        ORDER BY A.significant DESC, A.l2fold ASC
        '''

        return self.getAll(statement)


class SleuthResultsSig(SleuthResults):

    pattern = "(.*)_DEresults$"
    direction = ""
    where = "WHERE p_value NOT NULL AND significant == 1 AND ABS(l2fold) > 1"
