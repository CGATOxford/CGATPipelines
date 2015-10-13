from CGATReport.Tracker import TrackerSQL, SingleTableTrackerRows,\
    SingleTableTrackerHistogram


class ProjectTracker(TrackerSQL):
    """Base class for trackers in this report"""


class PicardDuplicatesMetrics(ProjectTracker, SingleTableTrackerRows):
    table = "picard_duplicates_duplicate_metrics"


class PicardDuplicatesHistogram(ProjectTracker, SingleTableTrackerHistogram):
    table = "picard_duplicates_duplicate_histogram"
    column = "duplicates"


class TagCountsSummary(ProjectTracker, SingleTableTrackerRows):
    table = "counts_stats"
    fields = ("metric", )


class WindowsSummary(ProjectTracker, SingleTableTrackerRows):
    table = "windows_stats"
    fields = ("data", )


class WindowsSizes(ProjectTracker, SingleTableTrackerHistogram):
    table = "windows_hist"
    column = "residues"


class TagCountsCorrelations(ProjectTracker):
    pattern = "(.*)_correlation"

    def __call__(self, track):
        return self.getDict("SELECT * FROM %(track)s_correlation")


class TagCountsDesignSummary(ProjectTracker):
    pattern = "design(.*)_stats"

    def __call__(self, track):
        return self.getAll("SELECT * FROM design%(track)s_stats")
