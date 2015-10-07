import os
import glob
from collections import OrderedDict as odict
from CGATReport.Tracker import TrackerSQL, SingleTableTrackerRows
from CGATReport.Utils import PARAMS as P
from CGATReport.Utils import ResultBlock, ResultBlocks, layoutBlocks
import CGATPipelines.PipelineTracks as PipelineTracks

# parameterization
EXPORTDIR = P.get('readqc_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('readqc_datadir', P.get('datadir', '.'))
DATABASE = P.get('readqc_backend', P.get('sql_backend', 'sqlite:///./csvdb'))
PROCESSEDDIR = "processed.dir"

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.sra" % DATADIR), "(\S+).sra") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % DATADIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % DATADIR), "(\S+).fastq.1.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("*.csfasta.gz"), "(\S+).csfasta.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % PROCESSEDDIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % PROCESSEDDIR), "(\S+).fastq.1.gz")


UNPROCESSED_TRACKS = [x for x in TRACKS if not str(x).startswith("trimmed")]


class ReadqcTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class TrackerFastQC(ReadqcTracker):
    tracks = ["all"]

    def __call__(self, track, slice=None):

        edir = EXPORTDIR

        toc_text = []
        link_text = []

        tracks = sorted([x.asFile() for x in TRACKS])
        for track in tracks:

            for x, fn in enumerate(glob.glob(os.path.join(
                    EXPORTDIR, "fastqc", "%s*_fastqc" % track))):
                y = x + 1
                toc_text.append("* %(track)s-%(y)i_" % locals())
                link_text.append(
                    ".. _%(track)s-%(y)i: %(fn)s/fastqc_report.html" %
                    locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)
        rst_text = "\n%(toc_text)s\n\n%(link_text)s\n" % locals()
        return odict((("text", rst_text),))


class FilteringSummary(SingleTableTrackerRows):
    table = "filtering_summary"


class FastQCDetails(ReadqcTracker):
    tracks = ["all"]
    slices = ("adapter_content",
              "duplication_levels",
              "kmer_profiles",
              "per_base_n_content",
              "per_base_quality",
              "per_base_sequence_content",
              "per_sequence_gc_content",
              "per_sequence_quality",
              "per_tile_quality",
              "sequence_length_distribution")

    def __call__(self, track, slice=None):

        # note there are spaces behind the %(image)s directive to accommodate
        # for path substitution
        block = """.. figure:: %(image)s                                                                   
   :height: 300
"""
        result = odict()

        for track in sorted([x.asFile() for x in TRACKS]):
            if track.startswith("trimmed"):
                continue

            files = glob.glob(
                os.path.join(EXPORTDIR, "fastqc", "*%s*_fastqc" % track))

            text = ["%s\n\n" % track]

            blocks = ResultBlocks()

            for fn in sorted(files):

                image = os.path.abspath(
                    os.path.join(fn, "Images", "%s.png" % slice))
                if not os.path.exists(image):
                    continue

                blocks.append(ResultBlock(block % locals(),
                                          os.path.basename(fn)))

            text.append("\n".join(layoutBlocks(
                blocks, "grid")))

            result[track] = {'rst': "\n".join(text)}

        return result


class FastqcSummaryFull(ReadqcTracker):
    pattern = "(.*)_Basic_Statistics"
    slices = ("File type", "Filename", "Encoding",
              "Total Sequences", "Sequence Length", "%GC")

    def __call__(self, track, slice):
        return self.getAll(
            """SELECT * FROM %(track)s_Basic_Statistics
            WHERE measure = '%(slice)s'""")


class FastqcSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "basic_statistics_summary"


class OverRepresentedSequences(ReadqcTracker):
    pattern = "(.*)_Overrepresented_sequences"
    
    def __call__(self, track):
        return self.getAll(
            "SELECT * FROM %(track)s_Overrepresented_sequences")


class ProcessingComparison(ReadqcTracker):
    """compare before after for read-processing."""
    table = "basic_statistics_summary"
    tracks = sorted([x.asFile() for x in UNPROCESSED_TRACKS])
    slices = ("Total_Sequences", "Sequence_length")

    def __call__(self, track, slice):

        untrimmed = self.getValue(
            "SELECT %(slice)s FROM %(table)s WHERE track = '%(track)s'")

        trimmed = self.getValue(
            "SELECT %(slice)s FROM %(table)s "
            "WHERE track = 'trimmed-%(track)s'")

        return odict((("untrimmed", untrimmed),
                      ("trimmed", trimmed)))


class ProcessingDetails(ReadqcTracker):

    '''return summary of the read processing steps.'''
    pattern = "(.*)_processed$"

    def __call__(self, track):
        return self.getAll(
            """SELECT pair,input,output,pair,
            100.0 * output / input as percent
            FROM %(track)s_processed""")


class ProcessingSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "processing_summary"


