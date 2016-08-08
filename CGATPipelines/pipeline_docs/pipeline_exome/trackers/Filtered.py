import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class recessives(ExomeTracker):

    pattern = "(.*)_recessive_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT r.CHROM, r.POS, r.REF, r.ALT, r.ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, EFF, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, AC_Adj, AN_Adj, (CAST(AC_Adj AS FLOAT)/AN_Adj) AS ExAC, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, CLNDBN, CLNORIGIN, CLNSIG, FILTER, BaseQRankSum, FS, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_recessive_table AS r INNER JOIN all_samples_snpeff_table AS s ON r.CHROM = s.CHROM AND r.POS = s.POS WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND (ExAC is null OR ExAC<0.01) AND FILTER='PASS' " % locals())


class dominants(ExomeTracker):

    pattern = "(.*)_dominant_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT d.CHROM, d.POS, d.REF, d.ALT, d.ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, EFF, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, AC_Adj, AN_Adj, (CAST(AC_Adj AS FLOAT)/AN_Adj) AS ExAC, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, CLNDBN, CLNORIGIN, CLNSIG, FILTER, BaseQRankSum, FS, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_dominant_table AS d INNER JOIN all_samples_snpeff_table AS s ON d.CHROM = s.CHROM AND d.POS = s.POS WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND (ExAC<0.001 OR ExAC is null) AND FILTER='PASS' " % locals())


class deNovos(ExomeTracker):

    pattern = "(.*)_filtered_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT f.CHROM, f.POS, f.REF, f.ALT, f.ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, EFF, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, AC_Adj, AN_Adj, (CAST(AC_Adj AS FLOAT)/AN_Adj) AS ExAC, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, CLNDBN, CLNORIGIN, CLNSIG, FILTER, BaseQRankSum, FS, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_filtered_table AS f INNER JOIN all_samples_snpeff_table AS s ON f.CHROM = s.CHROM AND f.POS = s.POS WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND (ExAC<0.001 OR ExAC is null)  AND FILTER='PASS' " % locals())


class comp_hets(ExomeTracker):

    pattern = "(.*)_compound_hets_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT c.rowid, g.CHROM, g.POS, g.REF, g.ALT, g.ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, AC_Adj, AN_Adj, (CAST(AC_Adj AS FLOAT)/AN_Adj) AS ExAC, dbNSFP_ESP6500_EA_AF, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, g.CLNDBN, g.CLNORIGIN, g.CLNSIG, g.FILTER, s.EFF, g.BaseQRankSum, g.FS, g.MQ, g.MQ0, g.MQRankSum, g.QD, g.ReadPosRankSum, gene, family_genotypes, family_members FROM %(track)s_compound_hets_table AS c INNER JOIN all_samples_snpsift_table AS g ON c.gene = g.SNPEFF_GENE_NAME AND c.qual = g.QUAL AND c.depth = g.DP AND c.ref = g.REF AND c.alt = g.ALT INNER JOIN all_samples_snpeff_table AS s ON g.CHROM = s.CHROM AND g.POS = s.POS AND g.REF = s.REF AND g.ALT = s.ALT WHERE c.rowid IN (SELECT MAX(%(track)s_compound_hets_table.rowid) FROM %(track)s_compound_hets_table GROUP BY chrom, start, end, codon_change) AND (g.FILTER = 'PASS' OR g.FILTER = 'GENE_OF_INTEREST')" % locals())
