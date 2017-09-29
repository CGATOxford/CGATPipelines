from cpgReport import *
from CGATReport.Tracker import *

##########################################################################


class replicatedIntervalTranscriptOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding transcripts """
    mPattern = "_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_transcript_extended_pover1"

##########################################################################


class replicatedIntervalGeneOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding genes """
    mPattern = "_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_gene_extended_pover1"
