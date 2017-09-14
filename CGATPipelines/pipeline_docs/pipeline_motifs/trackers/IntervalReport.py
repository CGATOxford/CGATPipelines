
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

EXPORTDIR = P.get('motifs_exportdir',
                  P.get('report_exportdir', "export"))
DATADIR = P.get('motifs_datadir',
                P.get('report_datadir', "."))
DATABASE = P.get('motifs_backend',
                 P.get('report_backend', "sqlite:///./csvdb"))


class IntervalTracker(TrackerSQL):
    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self,
                            *args,
                            backend=DATABASE,
                            **kwargs)
