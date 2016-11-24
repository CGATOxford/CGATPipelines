from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram
from MappingReport import MappingTracker


class MappingSummary(MappingTracker, SingleTableTrackerRows):
    table = "view_mapping"


class PairedMappingSummary(MappingTracker, SingleTableTrackerRows):
    table = "view_mapping"
    # select only tracks with paired reads
    where = "pairs_total > 0"


class TophatSummary(MappingTracker, SingleTableTrackerRows):
    table = "tophat_stats"


class StarSummary(MappingTracker, SingleTableTrackerRows):
    table = "star_stats"


class BamSummary(MappingTracker, SingleTableTrackerRows):
    table = "bam_stats"


class PicardSummary(MappingTracker, SingleTableTrackerRows):
    table = "picard_stats_alignment_summary_metrics"


class PicardDuplicationSummary(MappingTracker, SingleTableTrackerRows):
    table = "picard_duplication_metrics"


class PicardRnaMetrics(MappingTracker, SingleTableTrackerRows):
    table = "picard_rna_metrics"


class PicardAlignmentSummaryMetrics(MappingTracker, SingleTableTrackerRows):
    table = "picard_stats_alignment_summary_metrics"


class PicardInsertSizeMetrics(MappingTracker, SingleTableTrackerRows):
    table = "picard_stats_insert_size_metrics"


class PicardDuplicatesMetrics(MappingTracker, SingleTableTrackerRows):
    table = "picard_duplicates_duplicate_metrics"


class PicardInsertSizeHistogram(MappingTracker, SingleTableTrackerHistogram):
    table = "picard_stats_insert_size_histogram"
    column = "insert_size"


class PicardDuplicatesHistogram(MappingTracker,
                                SingleTableTrackerHistogram):
    table = "picard_duplicates_duplicate_histogram"
    column = "duplicates"


class PicardQualityByCycleHistogram(MappingTracker,
                                    SingleTableTrackerHistogram):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"


class PicardQualityDistributionHistogram(MappingTracker,
                                         SingleTableTrackerHistogram):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"


class DuplicationMetricsTable(MappingTracker, SingleTableTrackerHistogram):

    table = "picard_complexity_histogram"

    def __call__(self, track=None, slice=None):
        cols = self.getColumns(self.table)
        if len(cols) == 0:
            return None

        fields = ", ".join(cols)
        data = self.getAll(
            "SELECT %(fields)s FROM %(table)s ORDER BY coverage_multiple")
        return data


class RnaBiasTable(MappingTracker, SingleTableTrackerHistogram):

    table = "picard_rna_histogram"
    column = "coverage_multiple"


class MappingFlagsMismatches(MappingTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nm"
    column = "nm"


class MappingFlagsHits(MappingTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nh"
    column = "nh"


class MappingQuality(MappingTracker, SingleTableTrackerHistogram):
    table = "bam_stats_mapq"
    column = "mapq"


class AlignmentQualityByCycle(MappingTracker, SingleTableTrackerHistogram):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"


class DuplicationMetrics(MappingTracker, SingleTableTrackerHistogram):
    table = "picard_complexity_histogram"
    column = "coverage_multiple"


class AlignmentQualityDistribution(MappingTracker,
                                   SingleTableTrackerHistogram):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"


class MappingContext(MappingTracker, SingleTableTrackerRows):
    table = "context_stats"


class FilteringSummary(MappingTracker, SingleTableTrackerRows):
    table = "mapping_stats"


class BigwigSummary(MappingTracker, SingleTableTrackerRows):
    table = "bigwig_stats"
