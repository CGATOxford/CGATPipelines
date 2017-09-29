from CGATReport.Tracker import *
import collections


# classes for returning data and plots for genome coverage analysis

class GenomeCoverage(TrackerSQL):

    def __call__(self, track, slice=None):

        result = collections.defaultdict(list)
        cols = ["observed", "expected"]
        for col in cols:
            for data in self.execute("""SELECT %(col)s FROM genome_coverage_observed_expected""" % locals()):
                result[col].append(data[0])
        return result
