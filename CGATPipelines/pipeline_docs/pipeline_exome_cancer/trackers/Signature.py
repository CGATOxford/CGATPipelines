from CGATReport.Tracker import *
from exomeReport import *


class SnpSummary(ExomeTracker):

    def __call__(self, track, slice=None):
        table = "mutational_signature"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        return self.getAll(statement)


class SnpSummaryTable(ExomeTracker):

    def __call__(self, track, slice=None):
        table = "mutational_signature_table"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        return self.getAll(statement)
