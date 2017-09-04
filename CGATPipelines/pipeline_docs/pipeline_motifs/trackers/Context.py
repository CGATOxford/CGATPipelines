
from CGATReport.Tracker import *
from IntervalReport import *


class ContextSummary(IntervalTracker, SingleTableTrackerRows):
    table = "context_stats"
