from CGATReport.Tracker import Status
import glob
import re


class ComparisonStatus(Status):
    '''pipeline status
    '''

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM md5_compare")
        return tuple([x[0] for x in d])

    slices = ('FileComparison',)

    def testFileComparison(self, track):
        '''
        PASS: All files exist and are the same.

        FAIL: There are extra files or missing files

        WARN: All files exist, but there are differences in the files.

        The value indicates the number of files missing, exta
        or different.
        '''

        data = dict(self.getRow("""SELECT missing, extra, different
        FROM md5_compare WHERE track = '%(track)s'"""))

        total = sum(data.values())

        if total == 0:
            status = "PASS"
        elif data['missing'] == 0 and data['extra'] == 0:
            status = "WARN"
        else:
            status = "FAIL"

        value = ",".join(["%s:%i" % (x, y) for x, y in data.items()])

        return status, value


class PipelineStatus(Status):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def testPipelineCompletion(self, track):
        '''Check if the pipeline started and completed successfully.

        PASS: the pipeline started and completed successfully.

        FAIL: the pipeline started, but did not complete successfully

        NA: the pipeline was not started.
        '''

        try:
            lines = open(track + ".log").readlines()
        except IOError:
            return 'NA', 'no log file'

        started = "not started"

        if len(lines) < 1:
            return 'FAIL', started

        x = re.search("# job started at ([^-]*) on", lines[1])
        if x:
            started = x.groups()[0]

        x = re.search(
            "# job finished in (\d+) seconds at ([^-]*) -- ", lines[-1])
        if not x:
            return 'FAIL', started
        else:
            return 'PASS', x.groups()[1]

    def testReportCompletion(self, track):
        '''Check if building the pipeline report completed successfully.

        PASS: the pipeline report was build successfully.

        FAIL: the pipeline report was not built successfully.

        NA: the pipeline report was not built.
        '''

        try:
            lines = open(track + ".report").readlines()
        except IOError:
            return 'NA', 'no report'

        lines = open(track + ".report").readlines()
        started = [x for x in lines if x.startswith("# job started")]
        finished = [x for x in lines if x.startswith("# job finished")]
        error = [x for x in lines if "ERROR" in lines]

        if len(started) == 0:
            return 'WARN', 'never started'

        if len(finished) == 0:
            return 'FAIL', 'started, but never finished'

        if error:
            return 'FAIL', 'report caused errors'
        else:
            return 'PASS', 'report completed'
