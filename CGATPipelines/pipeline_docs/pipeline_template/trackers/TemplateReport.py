
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P


###################################################################
###################################################################
# parameterization

EXPORTDIR = P['template_exportdir']
DATADIR = P['template_datadir']
DATABASE = P['template_backend']

###########################################################################


class TemplateTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
