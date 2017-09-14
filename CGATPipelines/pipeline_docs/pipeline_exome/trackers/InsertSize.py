
from CGATReport.Tracker import *
from exomeReport import *


class PicardInsertSizeStats(ExomeTracker, SingleTableTrackerRows):
    table = "picard_isize_stats"
