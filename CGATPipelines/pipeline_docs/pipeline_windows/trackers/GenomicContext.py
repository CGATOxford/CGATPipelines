from MedipReport import ProjectTracker
from CGATReport.Tracker import SingleTableTrackerRows, \
    SingleTableTrackerEdgeList


class ContextSummary(ProjectTracker, SingleTableTrackerRows):
    table = "context_stats"


class GatFold(ProjectTracker, SingleTableTrackerEdgeList):

    '''fold change matrix.'''
    row = "sample"
    column = "track"
    value = "fold"
    where = "pvalue < 0.05"


class GatLogFold(ProjectTracker):

    '''logfold - colour is signficance'''
    fdr = 0.05
    as_tables = True

    def __call__(self, track):
        return self.getDict("""SELECT track, l2fold,
        CASE WHEN qvalue < %(fdr)f THEN 'red' ELSE 'gray' END AS colour
        FROM %(track)s ORDER BY fold""")


class GatResults(ProjectTracker):
    as_tables = True

    def __call__(self, track):
        return self.getAll("SELECT * FROM %(track)s")


class GatSignificantResults(ProjectTracker):
    as_tables = True
    fdr = 0.05

    def __call__(self, track):
        return self.getAll("SELECT * FROM %(track)s WHERE qvalue < %(fdr)f")


class GatTableContext:
    pattern = "(.*)_context_gat$"


_gat_analysis = {"Results": GatResults,
                 "SignificantResults": GatSignificantResults,
                 "Fold": GatLogFold,
                 "LogFold": GatLogFold}

_gat_sets = {"Context": GatTableContext}

for a, aa in list(_gat_analysis.items()):
    for b, bb in list(_gat_sets.items()):
        n = "Gat%s%s" % (a, b)
        globals()[n] = type(n, (bb, aa), {})

