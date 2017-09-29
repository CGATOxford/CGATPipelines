from CGATReport.Tracker import TrackerSQL


class ProjectTracker(TrackerSQL):
    '''Define convenience tracks for plots'''


class WordFrequencies(ProjectTracker):
    pattern = "(.*)_counts"

    def __call__(self, track):
        return self.getValues("SELECT freq FROM %(track)s_counts")
