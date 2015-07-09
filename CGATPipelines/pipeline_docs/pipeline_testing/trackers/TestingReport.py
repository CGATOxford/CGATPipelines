import os
import glob
import collections
from collections import OrderedDict as odict

from CGATReport.Tracker import TrackerSQL
from CGATReport.Utils import PARAMS as P

###################################################################
###################################################################
# parameterization

EXPORTDIR = P.get('testing_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('testing_datadir', P.get('datadir', '.'))
DATABASE = P.get('testing_backend', P.get('sql_backend', 'sqlite:///./csvdb'))


class TestingTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

LogFileSummary = collections.namedtuple(
    "LogFileSummary",
    "report_warning report_error info debug warning error")


def summarizeLogFile(filename):
    '''summarize a CGATReport logfile.'''

    info, debug, warning, error = 0, 0, 0, 0
    report_error, report_warning = 0, 0

    with open(filename) as f:
        for line in f:
            words = line.split()
            if len(words) < 3:
                continue
            w = words[2]
            if w == "INFO":
                info += 1
            elif w == "DEBUG":
                debug += 1
            elif w == "WARNING":
                if words[3].startswith("CGATReport-Warning"):
                    report_warning += 1
                else:
                    warning += 1
            elif w == "ERROR":
                if words[3].startswith("CGATReport-Error"):
                    report_error += 1
                else:
                    error += 1

    return LogFileSummary._make((report_warning, report_error,
                                 info, debug,
                                 warning, error))


class ReportTable(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def __call__(self, track):

        try:
            logfileresult = summarizeLogFile(
                os.path.join(track + ".dir", "cgatreport.log"))
        except IOError:
            return

        report_file = os.path.join(track + ".dir", "report.log")
        fn = os.path.abspath(
            os.path.join(track + ".dir", "report", "html", "contents.html"))

        r = odict()
        r["link"] = "`%(track)s <%(fn)s>`_" % locals()
        r["report_error"] = logfileresult.report_error
        r["report_warning"] = logfileresult.report_warning
        r["error"] = logfileresult.error
        r["warning"] = logfileresult.warning
        return r


class XReportTable(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def __call__(self, track):

        logfileresult = summarizeLogFile(
            os.path.join(track + ".dir", "cgatreport.log"))
        report_file = os.path.join(track + ".dir", "report.log")

        toc_text = []
        link_text = []

        fn = os.path.join(track + ".dir", "report", "html", "contents.html")
        toc_text.append("* %(track)s_" % locals())
        link_text.append(".. _%(track)s: %(fn)s" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))


class FilesWithProblems(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    slices = ("files_different_lines",
              "files_different_md5",
              "files_extra", "files_missing")

    def __call__(self, track, slice):

        statement = """SELECT %(slice)s FROM md5_compare
        WHERE track = '%(track)s'""" % locals()

        data = self.getValue(statement)
        if data is None:
            return

        files = data.split(",")
        # do not use :download: as
        # that will include the file in the report.
        # and thus wastes space.
        return ['`%s <file://%s>`_' %
                (f,
                 os.path.abspath(f)) for f in files]
