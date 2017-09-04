
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P



###################################################################
###################################################################
# parameterization
EXPORTDIR = P.get('intervals_exportdir',
                  P.get('exportdir', 'export'))
DATADIR = P.get('intervals_datadir',
                P.get('datadir', '.'))
DATABASE = P.get('intervals_backend',
                 P.get('sql_backend', 'sqlite:///./csvdb'))


class IntervalTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
