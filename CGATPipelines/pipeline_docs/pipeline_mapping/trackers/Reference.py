import math
from collections import OrderedDict as odict
from MappingReport import MappingTracker


class ReferenceData(MappingTracker):

    """Base class f or Trackers accessing reference table."""
    pattern = "(.*)_transcript_counts$"
    reference = "refcoding"


class TranscriptCoverage(ReferenceData):

    """Coverage of reference transcripts."""
    mXLabel = "overlap / %"

    def __call__(self, track, slice=None):
        return self.getValues(
            """SELECT coverage_sense_pcovered
            FROM %(track)s_transcript_counts
            WHERE coverage_sense_nval > 0""")


class GeneCoverage(ReferenceData):

    '''Coverage of reference genes - max transcript coverage per gene.'''

    def __call__(self, track, slice=None):
        return self.getValues(
            """SELECT max(c.coverage_sense_pcovered) FROM
            %(track)s_transcript_counts as c,
            %(reference)s_transcript2gene as i
            WHERE c.coverage_sense_nval > 0
            AND i.transcript_id = c.transcript_id
            GROUP BY i.gene_id""")


class MeanVsMaxReadDepth(ReferenceData):

    """maxmimum read depth versus mean read depth of :term:`reference`
    genes.  Dots are coloured by the log(length) of a
    :term:`reference` gene.

    """

    mXLabel = "mean read depth"
    mYLabel = "maximum read depth"

    def __call__(self, track, slice=None):
        reference = self.reference
        statement = """
        SELECT coverage_sense_mean, coverage_sense_max, exons_sum
        FROM %(track)s_transcript_counts""" % locals()

        data = [(x[0], x[1], math.log(x[2]))
                for x in self.get(statement) if x[2] > 0]
        return odict(list(zip(("mean coverage", "max coverage", "length"),
                              list(zip(*data)))))


class MeanVsMedianReadDepth(ReferenceData):

    """maxmimum read depth versus mean read depth of :term:`reference`
    genes.  Dots are coloured by the log(length) of a
    :term:`reference` gene.

    """

    mXLabel = "mean read depth"
    mYLabel = "median read depth"

    def __call__(self, track, slice=None):
        reference = self.reference
        statement = """
        SELECT coverage_sense_mean, coverage_sense_median, exons_sum
        FROM %(track)s_transcript_counts""" % locals()

        data = [(x[0], x[1], math.log(x[2]))
                for x in self.get(statement) if x[2] > 0]
        return odict(list(zip(("mean coverage", "median coverage", "length"),
                              list(zip(*data)))))


class ReadDirectionality(MappingTracker):

    '''return antisense / sense direction of reads in introns/genes.

    +1 is added as pseudo-count.
    '''

    pattern = "(.*)_intron_counts$"

    slices = ("intron", "transcript")

    def __call__(self, track, slice=None):
        data = self.getValues(
            """SELECT CAST( (coverage_antisense_nreads + 1) AS FLOAT) /
            (coverage_sense_nreads + 1)
            FROM %(track)s_%(slice)s_counts """)
        return odict((("direction", data),))


class IntronicExonicReadDepth(MappingTracker):

    '''return the maximum read depth in introns
    and exons of a gene.

    +1 is added as pseudo-count.
    '''
    pattern = "(.*)_intron_counts$"

    slices = ("anysense", "antisense", "sense")
    min_coverage = 0

    def __call__(self, track, slice=None):
        data = self.getAll(
            """SELECT e.coverage_%(slice)s_max + 1 AS exon,
            i.coverage_%(slice)s_max + 1 as intron
            FROM %(track)s_gene_counts as e, %(track)s_intron_counts as i
            WHERE e.gene_id = i.gene_id
            AND e.coverage_%(slice)s_max >= %(min_coverage)i""")
        return data
