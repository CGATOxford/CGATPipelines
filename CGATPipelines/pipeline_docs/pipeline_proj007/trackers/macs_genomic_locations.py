import cpgReport


##########################################################################


class genomicFeatures(cpgReport.cpgTracker):

    """return overlap of interval with genomic features """

    mPattern = "_merged_genomic_features$"

    def __call__(self, track, slice=None):
        data = self.getAll( """SELECT feature_class, count(distinct gene_id) as intervals FROM (
                               SELECT gene_id,
                               CASE WHEN  tss_extended_pover1 > 0  THEN 'TSS'
                               WHEN genes_pover1 > 0 THEN 'Gene'
                               WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                               WHEN downstream_flank_pover1 THEN 'Downstream'
                               ELSE 'Intergenic'
                               END AS feature_class
                               FROM %(track)s_merged_genomic_features)
                               group by feature_class""" % locals() )
        return data
