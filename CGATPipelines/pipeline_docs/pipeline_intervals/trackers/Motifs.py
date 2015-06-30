import os
import re
import math
import numpy
import numpy.ma
import xml.etree.ElementTree
from IntervalReport import *

##########################################################################
##########################################################################
##########################################################################
# compute MAST curve
##########################################################################


def computeMastCurve(evalues):
    '''compute a MAST curve.

    see http://www.nature.com/nbt/journal/v26/n12/extref/nbt.1508-S1.pdf

    returns a tuple of arrays (evalues, with_motifs, explained )
    '''

    if len(evalues) == 0:
        raise ValueError("no data")

    mi, ma = math.floor(min(evalues)), math.ceil(max(evalues))

    if mi == ma:
        raise ValueError("not enough data")

    hist, bin_edges = numpy.histogram(evalues, bins=numpy.arange(mi, ma, 1.0))
    with_motifs = numpy.cumsum(hist)
    explained = numpy.array(with_motifs)

    for x, evalue in enumerate(bin_edges[:-1]):
        explained[x] -= evalue

    explained[explained < 0] = 0

    return bin_edges[:-1], with_motifs, explained

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


def getFDR(samples, control, num_bins=1000):
    '''return the score cutoff at a certain FDR threshold using scores and
    control scores.  Note that this method assumes that a higher score
    is a better result.

    The FDR is defined as fdr = expected number of false positives
    (FP) / number of positives (P)

    Given a certain score threshold s , the following will be used as
    approximations:

    FP: the number of controls with a score of less than or equal to s. These
        are all assumed to be false positives. Both samples and control should 
        contain rougly equal number of entries, but FP is scaled to be equivalent to P. 

    P: the number of samples with a score of less than or equal to s. These
       are a mixture of both true and false positives.

    returns the score cutoff at FDR threshold.

    '''

    if len(samples) == 0 or len(control) == 0:
        return None, None

    mi1, ma1 = min(samples), max(samples)
    mi2, ma2 = min(control), max(control)
    mi, ma = min(mi1, mi2), max(ma1, ma2)
    hist_samples, bin_edges_samples = numpy.histogram(
        samples, range=(mi, ma), bins=num_bins)
    hist_control, bin_edges_control = numpy.histogram(
        control, range=(mi, ma), bins=num_bins)
    hist_samples = hist_samples[::-1].cumsum()
    hist_control = hist_control[::-1].cumsum()
    bin_edges = bin_edges_samples[::-1]

    # correct for unequal size in the two sets
    correction = float(len(samples)) / len(control)

    fdrs = []
    m = 0
    for s, p, fp in zip(bin_edges[:-1], hist_samples, hist_control):
        if p != 0:
            fdr = min(1.0, correction * float(fp) / p)
            m = max(m, fdr)
        fdrs.append(m)

    return bin_edges[:-1][::-1], fdrs[::-1]

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


def getMastEvalueCutoff(evalues, control_evalues, fdr=0.1):
    '''return the E-Value cutoff at a certain FDR threshold
    using the control tracks.

    returns the evalue cutoff at FDR threshold.
    '''
    bin_edges, fdrs = getFDR(
        [-x for x in evalues], [-x for x in control_evalues])

    for bin, f in zip(bin_edges, fdrs):
        if f < fdr:
            return -bin

    return 0


##########################################################################
##########################################################################
##########################################################################
# Base class for mast analysis
##########################################################################
class Mast(IntervalTracker):
    pattern = "(.*)_mast$"

    def getSlices(self):
        return self.getValues("SELECT DISTINCT motif FROM motif_info")

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class MastFDR(Mast):

    """return arrays of glam2scan scores
    """
    pattern = "(*.)_mast$"

    def __call__(self, track, slice=None):
        evalues = self.getValues(
            "SELECT -evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        control_evalues = self.getValues(
            "SELECT -min_evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        bin_edges, fdrs = getFDR(evalues, control_evalues)
        if bin_edges is None:
            return odict()
        bin_edges = [-x for x in bin_edges]
        print len(bin_edges), len(fdrs)

        return odict((("score", bin_edges),
                      ("fdr", fdrs)))

##########################################################################
##########################################################################
##########################################################################
# Annotation of bases with SNPs
##########################################################################


class MastSummary(Mast):

    """return summary of mast results.

    Return for each track the number of intervals in total,
    the number of intervals submitted to mast, 

    The evalue used as a MAST curve cutoff, the number and % explained using the Mast cutoff.

    """

    mEvalueCutoff = 1
    mFDR = 0.1

    def __call__(self, track, slice=None):

        data = []
        nintervals = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_intervals" % locals())
        data.append(("nintervals", nintervals))
        data.append(("nmast", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())))

        evalues = self.getValues(
            "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        if len(evalues) <= 1:
            return odict()

        try:
            bin_edges, with_motifs, explained = computeMastCurve(evalues)
        except ValueError, msg:
            return odict((("msg", msg),))

        if len(explained) == 0:
            return odict((("msg", "no data"), ))

        am = numpy.argmax(explained)
        evalue = bin_edges[am]

        intervals_and_peakvals = \
            self.get("""
            SELECT peakval, evalue
            FROM %(track)s_mast AS m, %(track)s_intervals AS i
            WHERE i.interval_id = m.id AND motif = '%(slice)s'
            ORDER BY peakval""" % locals())

        intervals_with_motifs = len(
            [x for x in intervals_and_peakvals if x[1] <= evalue])
        ntop = nintervals / 4
        top25_with_motifs = len(
            [x for x in intervals_and_peakvals[-ntop:] if x[1] <= evalue])
        bottom25_with_motifs = len(
            [x for x in intervals_and_peakvals[:ntop] if x[1] <= evalue])

        data.append(("MC-Evalue", bin_edges[am]))
        data.append(("MC-explained", explained[am]))
        data.append(
            ("MC-explained / %", "%5.2f" % (100.0 * explained[am] / nintervals)))
        data.append(("MC-with-motif", intervals_with_motifs))
        data.append(("MC-with-motif / %", "%5.2f" %
                    (100.0 * intervals_with_motifs / nintervals)))
        data.append(("MC-top25-with-motif", top25_with_motifs))
        if ntop == 0:
            ntop = 1
        data.append(("MC-top25-with-motif / %", "%5.2f" %
                    (100.0 * top25_with_motifs / ntop)))
        data.append(("MC-bottom25-with-motif", bottom25_with_motifs))
        data.append(("MC-bottom25-with-motif / %", "%5.2f" %
                    (100.0 * bottom25_with_motifs / ntop)))

        # use control intervals to compute FDR
        control_evalues = self.getValues(
            "SELECT min_evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        evalue_cutoff = getMastEvalueCutoff(
            evalues, control_evalues, self.mFDR)
        nexplained = len([x for x in evalues if x <= evalue_cutoff])
        data.append(("FDR-Evalue", evalue_cutoff))
        data.append(("FDR-explained", nexplained))
        data.append(
            ("FDR-explained / %", "%5.2f" % (100.0 * nexplained / nintervals)))

        # use a pre-defined E-value threshold
        threshold = self.mEvalueCutoff
        data.append(("Evalue", self.mEvalueCutoff))
        n = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s' AND evalue <= %(threshold)f" % locals())
        data.append(("threshold-explained", n))
        data.append(
            ("threshold-explained / %", "%5.2f" % (100.0 * n / nintervals)))

        # use no threshold (nmatches > 1)
        n = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s' AND nmatches > 0" % locals())
        data.append(("nmatches-explained", n))
        data.append(
            ("nmatches-explained / %", "%5.2f" % (100.0 * n / nintervals)))

        intervals_and_peakvals = \
            self.get("""
                               SELECT peakval, nmatches
                               FROM %(track)s_mast AS m, %(track)s_intervals AS i
                               WHERE i.interval_id = m.id AND motif = '%(slice)s'
                               ORDER BY peakval""" % locals())

        intervals_with_motifs = len(
            [x for x in intervals_and_peakvals if x[1] > 0])
        ntop = nintervals / 4
        top25_with_motifs = len(
            [x for x in intervals_and_peakvals[-ntop:] if x[1] > 0])
        bottom25_with_motifs = len(
            [x for x in intervals_and_peakvals[:ntop] if x[1] > 0])

        if ntop == 0:
            ntop = 1
        data.append(("nmatches-top25-with-motif", top25_with_motifs))
        data.append(("nmatches-top25-with-motif / %", "%5.2f" %
                    (100.0 * top25_with_motifs / ntop)))
        data.append(("nmatches-bottom25-with-motif", bottom25_with_motifs))
        data.append(("nmatches-bottom25-with-motif / %", "%5.2f" %
                    (100.0 * bottom25_with_motifs / ntop)))

        return odict(data)


class MastQuickSummary(Mast):

    """return a quicker summary of MC analysis of mast results.

    Return for each track the number of intervals in total, the number
    of intervals submitted to mast, The evalue used as a MAST curve
    cutoff, the number and % explained using the Mast cutoff.

    """

    mEvalueCutoff = 1
    mFDR = 0.1

    def __call__(self, track, slice=None):

        data = []
        nintervals = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_intervals" % locals())
        data.append(("nintervals", nintervals))

        evalues = self.getValues(
            "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        if len(evalues) <= 1:
            return odict()

        try:
            bin_edges, with_motifs, explained = computeMastCurve(evalues)
        except ValueError, msg:
            return odict((("msg", msg),))

        if len(explained) == 0:
            return odict((("msg", "no data"), ))

        am = numpy.argmax(explained)
        evalue = bin_edges[am]

        data.append(("MC-With-Motifs", with_motifs[am]))
        data.append(("MC-Evalue", bin_edges[am]))
        data.append(("MC-explained", explained[am]))
        data.append(
            ("MC-explained / %", "%5.2f" % (100.0 * explained[am] / nintervals)))

        return odict(data)


class MastMotifEvalues(Mast):

    '''distribution of evalues.'''

    def __call__(self, track, slice=None):
        r = odict()
        for x in ("evalue", "l_evalue", "r_evalue"):
            r[x] = self.getValues(
                "SELECT %(x)s FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        return r


class MastMotifLocation(Mast):

    '''plot median position of motifs versus the peak location.

    The position is centered around 1kb.
    '''

    def __call__(self, track, slice=None):

        data = numpy.array(self.getValues(
            """SELECT (i.peakcenter -
            (m.start + (m.end - m.start) / 2)) / 500.0
            FROM %(track)s_mast as m,
                 %(track)s_intervals as i
            WHERE i.interval_id = m.id AND motif = '%(slice)s'
            AND m.nmatches = 1"""
            % locals()), numpy.float)

        data[data < -1.0] = -1.0
        data[data > 1.0] = 1.0

        return odict((("distance", data),))


class MastMotifLocationMiddle(Mast):

    '''plot median position of motifs versus the center of the interval.'''

    def __call__(self, track, slice=None):

        # difference between
        #   middle of interval: i.start + i.length / 2
        #   middle of motif: m.start + (m.end - m.start) / 2
        # divide by (intervalsize - motifsize) / 2
        #
        # only take single matches (multiple matches need not be centered)
        data = self.getValues( """SELECT ((i.start + i.length / 2) - (m.start + (m.end - m.start) / 2)) 
                                         / ((CAST(i.length AS FLOAT) - (m.end - m.start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.nmatches = 1"""
                               % locals())
        return odict((("distance", data),))


class MastControlLocationMiddle(Mast):

    '''plot median position of controls versus the center of the interval.'''

    def __call__(self, track, slice=None):

        data1 = self.getValues( """SELECT ( (m.r_length / 2) - (m.r_start + (m.r_end - m.r_start) / 2) ) / ((CAST( m.r_length as float) - (m.r_end - m.r_start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.r_nmatches = 1"""
                                % locals())
        data2 = self.getValues( """SELECT ( (m.l_length / 2) - (m.l_start + (m.l_end - m.l_start) / 2) ) / ((CAST( m.l_length as float) - (m.l_end - m.l_start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.l_nmatches = 1"""
                                % locals())

        return odict((("distance", data1 + data2),))


class MastCurve(Mast):

    """Summary stats of mast results.
    """
    pattern = "(.*)_mast$"

    def __call__(self, track, slice=None):

        evalues = self.getValues(
            "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())

        if len(evalues) == 0:
            return odict()
        try:
            bin_edges, with_motifs, explained = computeMastCurve(evalues)
        except ValueError, msg:
            return odict()

        data = odict()
        data["with_motifs"] = odict(
            (("evalue", bin_edges), ("with_motifs", with_motifs)))
        data["explained"] = odict(
            (("evalue", bin_edges), ("explained", explained)))

        return data




class MastPeakValWithMotif(Mast):

    '''return for each peakval the proportion of intervals
    that have a motif.

    This tracker uses the nmatches indicator as cutoff.
    '''

    def __call__(self, track, slice):

        # obtain evalue distribution
        data = self.get('''
        SELECT i.peakval, m.nmatches
        FROM %(track)s_intervals AS i,
             %(track)s_mast AS m
        WHERE m.id = i.interval_id \
           AND m.motif = '%(slice)s' ORDER BY i.peakval DESC''' % locals())

        result = Stats.getSensitivityRecall(
            [(int(x[0]), x[1] > 0) for x in data])

        return odict(zip(("peakval", "proportion with motif", "recall"),
                         zip(*result)))



class MemeInputSequenceComposition(IntervalTracker):

    '''distribution of sequence composition in sequences
       submitted to motif searches.'''
    pattern = "(.*)_motifseq_stats"
    slices = ('nA', 'nAT', 'nC', 'nG', 'nGC', 'nN', 'nT',
              'nUnk', 'pA', 'pAT', 'pC', 'pG', 'pGC', 'pN', 'pT')

    def __call__(self, track, slice):
        return self.getValues('''SELECT %(slice)s
                                 FROM %(track)s_motifseq_stats''')


class MotifRuns(IntervalTracker):

    def getSlices(self):
        methods = self.getValues("SELECT DISTINCT method FROM %s_summary" %
                              self.prog)
        return ["%s (:ref:`%s_details`)" % (method, method)
                for method in methods]

    def getTracks(self):
        slices = ",".join(self.getSlices())
        return self.getValues("SELECT DISTINCT track FROM %s_summary" %
                              self.prog)

    def getExportDir(self, track, slice):
        memedir = os.path.abspath(
            os.path.join(EXPORTDIR, slice + ".dir", "%s.%s" %
                         (track, self.prog)))
        return memedir

    def checkexists(self, track, slice):
        return os.path.exists(self.getExportDir(track, slice))

    def getNMotifs(self, track, slice):
        export_dir = self.getExportDir(track, slice)
        tree = xml.etree.ElementTree.ElementTree()
        tree.parse(os.path.join(export_dir, "%s.xml" % self.prog))
        motifs = tree.find("motifs")
        nmotifs = 0
        for motif in motifs.getiterator("motif"):
            nmotifs += 1
        return nmotifs

    def getTomTomLink(self, track, slice):
        tomtomdir = os.path.abspath(
            os.path.join(EXPORTDIR, "tomtom", "%s.tomtom" % track))
        if os.path.exists(tomtomdir):
            tomtomlink = "`tomtom_%s <%s/tomtom.html>`_" % (track, tomtomdir)
        else:
            tomtomlink = "NA"
        return tomtomlink

    def getInputInfo(self, track, slice):
        return []

    def getLink(self, track, slice):

        return "`%s_%s <%s/%s.html>`_" % (self.prog,
                                          track,
                                          self.getExportDir(track,slice),
                                          self.prog)

    def __call__(self, track, slice=None):

        slice = re.match("(.+) \(.+\)", slice).groups()[0]
        data = []
        if not self.checkexists(track, slice):
            return None

        data.extend(self.getInputInfo(track, slice))
        data.append(("nmotifs", self.getNMotifs(track, slice)))
        data.append(("Link",self.getLink(track, slice)))
        data.append(("TomTom",self.getTomTomLink(track, slice)))

        print data
        return odict(data)

       
class MemeRuns(MotifRuns):
    prog = "meme"
    
    def getInputInfo(self, track, slice):
        
        tree = xml.etree.ElementTree.ElementTree()
        tree.parse(os.path.join(self.getExportDir(track, slice), "meme.xml"))
        model = tree.find("model")
        return [("nsequences", int(model.find("num_sequences").text)),
                ("nbases", int(model.find("num_positions").text))]


class DremeRuns(MotifRuns):
    prog = "dreme"

    def getInputInfo(self, track, slice):

        track = re.sub("-", "_", track)

        try:
            positives, negatives = re.match("(.+)_vs_(.+)", track).groups()
        except:
            positives = negatives = track

        def _getinfo(track):
            n = self.getValue("SELECT COUNT (*) FROM %s_dreme_motifseq_stats"
                              % track)
            length = self.getValue(
                "SELECT SUM(nA+nG+nC+nT) FROM %s_dreme_motifseq_stats" % track)
            return (n, length)

        return zip(("n_positive_seqs", "n_positive_bases"), 
                   _getinfo(positives)) + \
            zip(("n_negative_seqs", "n_negatives_bases"),
                _getinfo(negatives))


class MemeChipRuns(MotifRuns):
    prog="memechip"

    def getInputInfo(self, track, slice):
        track = re.sub("-", "_", track)
        
        return [("nsequences", 
                 self.getValue(
                     "SELECT COUNT (*) from %s_memechip_motifseq_stats" % track)),
                ("nbases", 
                 self.getValue(
                     "SELECT SUM(nA+nG+nT+nC) FROM %s_memechip_motifseq_stats" % track))]

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM memechip_seeds")

    def getSlices(self):
        return ["MemeChip (:ref:`memechip_details`)"]

    def getNMotifs(self, track, slice=None):
        return self.getValue('''SELECT COUNT(*) from memechip_seeds 
                                WHERE track = '%(track)s' ''')

    def getLink(self, track, slice):
        return "`Full MemeChip Output <%s/index.html>`_" \
            % self.getExportDir(track, "memechip")

    def checkexists(self, track, slice):
        return True
    
    
def line2lines(values):
    output = []
    values = list(set(str(values).split()))
    
    if len(values) == 1:
        return values[0]
    while len(values) > 0:
        tmp = []
        while sum(len(x) for x in tmp) < 20 and len(values) > 0:
            tmp.append(values.pop(0))
            output.append(" ".join(tmp))
            
        output = "\n".join(output)
        
    return output


class MemeResults(IntervalTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM meme_summary")

    def getSlices(self):
        return self.getValues("SELECT DISTINCT method FROM meme_summary")

    def __call__(self, track, slice=None):

        statement = '''SELECT DISTINCT
                              primary_id as ID,
                              group_concat(track2, " ") as other_tracks
                      FROM %(slice)s_seeds as seeds
                      LEFT JOIN %(slice)s_track_comparisons as comps
                           ON seeds.track = comps.track1
                          AND seeds.primary_id = comps.Query_ID
                      WHERE track = '%(track)s'
                      GROUP BY primary_id'''

        alt_statement = '''SELECT DISTINCT primary_id as ID,
                                           'NA' as other_tracks
                                 FROM %(slice)s_seeds
                                 WHERE track = '%(track)s' '''
        try:
            motifs_to_keep = self.getDataFrame(statement)
        except:
            motifs_to_keep = self.getDataFrame(alt_statement)

        motifs_to_keep.set_index("ID", inplace=True)
        
        resultsdir = os.path.abspath(
            os.path.join(EXPORTDIR, slice + ".dir", "%s.meme" % track))
        if not os.path.exists(resultsdir):
            return None
 
        tree = xml.etree.ElementTree.ElementTree()
        tree.parse(os.path.join(resultsdir, "meme.xml"))

        motifs = tree.find("motifs")
        nmotif = 0
        result = odict()
        for motif in motifs.getiterator("motif"):
            nmotif += 1
            if nmotif not in motifs_to_keep.index:
                continue

            motif_img = "%s/logo%i.png" % (resultsdir, nmotif)
            motif_rc_img = "%s/logo_rc%i.png" % (resultsdir, nmotif)
            img, rc_img = "na", "na"
            if os.path.exists(motif_img):
                img = '''.. image:: %s
   :scale: 25%%''' % motif_img
            if os.path.exists(motif_rc_img):
                rc_img = '''.. image:: %s
   :scale: 25%%''' % motif_rc_img

            result[str(nmotif)] = odict((
                ("width", motif.get("width")),
                ("evalue", motif.get("e_value")),
                ("information content", motif.get("ic")),
                ("sites", motif.get("sites")),
                ("link", "`meme_%s_%i <%s/meme.html#summary%i>`_" %
                 (track, nmotif, resultsdir, nmotif)),
                ("img", img),
                ("rev", rc_img),
                ("Similar motif found in:", 
                 line2lines(motifs_to_keep.loc[nmotif]["other_tracks"])),
            ))

        return result


class DremeResults(IntervalTracker):
    
    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM dreme_summary")

    def getSlices(self):
        return self.getValues("SELECT DISTINCT method FROM dreme_summary")

    def __call__(self, track, slice=None):

        statement = '''SELECT DISTINCT
                              primary_id as ID,
                              group_concat(track2, " ") as other_tracks
                      FROM %(slice)s_seeds as seeds
                      LEFT JOIN %(slice)s_track_comparisons as comps
                           ON seeds.track = comps.track1
                          AND seeds.primary_id = comps.Query_ID
                      WHERE track = '%(track)s'
                      GROUP BY primary_id'''

        alt_statement = '''SELECT DISTINCT primary_id as ID,
                                           'NA' as other_tracks
                                 FROM %(slice)s_seeds
                                 WHERE track = '%(track)s' '''
        try:
            motifs_to_keep = self.getDataFrame(statement)
        except:
            print alt_statement % locals()
            motifs_to_keep = self.getDataFrame(alt_statement)

        motifs_to_keep.set_index("ID", inplace=True)

        resultsdir = os.path.abspath(
            os.path.join(EXPORTDIR, slice + ".dir", "%s.dreme" % track))
        if not os.path.exists(resultsdir):
            return None

        tree = xml.etree.ElementTree.ElementTree()
        tree.parse(os.path.join(resultsdir, "dreme.xml"))
        
        model = tree.find("model")
        num_positives = int(model.find("positives").get("count"))
        num_negatives = int(model.find("negatives").get("count"))

        result = odict()

        motifs = tree.find("motifs")
        nmotif = 0
        for motif in motifs.getiterator("motif"):

            seq = motif.get("seq")
            nmotif += 1
            id = motif.get("id")

            if seq not in motifs_to_keep.index.values:
                print "%s not in %s" % (seq, motifs_to_keep.index.values)
                continue

            motif_img = "%(resultsdir)s/%(id)snc_%(seq)s.png" % locals()
            rc_file = glob.glob(os.path.join(resultsdir, "%(id)src_*" % locals()))

            img, rc_img = "na", "na"
            if os.path.exists(motif_img):
                img = '''.. image:: %s
   :scale: 25%% ''' % motif_img
                
            if len(rc_file) > 0:
                rc_img = '''.. image:: %s
   :scale: 25%%''' % rc_file[0]

            p = float(motif.get("p"))
            n = float(motif.get("n"))
            
            try:
                enrichment = (p/num_positives)/(n/num_negatives)
                enrichment = "{:.0%}".format(
                    enrichment)
            except ZeroDivisionError:
                enrichment = "inf"

            result[str(nmotif)] = odict((
                ("sequence", seq),
                ("evalue", motif.get("evalue")),
                ("positives", "{:d}/{} ({:.0%})".format(int(p), num_positives,
                                                        p/num_positives)),
                ("negatives", "{:d}/{} ({:.0%})".format(int(n), num_negatives,
                                                        n/num_negatives)),
                ("enrichment", enrichment),
                ("link", "`dreme_%s <%s/dreme.html>`_" %
                 (track, resultsdir)),
                ("img", img),
                ("rc_img", rc_img),
                ("Simialr motif found in:", 
                 line2lines(motifs_to_keep.loc[seq]["other_tracks"])),
            ))

        return result
        
class MemeChipResults(IntervalTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM memechip_seeds")

    def __call__(self, track):

        print EXPORTDIR
        resultsdir = os.path.abspath(os.path.join(EXPORTDIR, "memechip.dir",
                                                  track + ".memechip"))

        statement = '''SELECT DISTINCT
                              primary_id as ID,
                              consensus,
                              nsites as '# seed matches',
                              nClustered as "# in cluster",
                              totalHits as '# cluster matches',
                              E as 'E-value',
                              group_concat(track2, " ") as "Similar motif found in"
                      FROM memechip_seeds as seeds
                      LEFT JOIN meme_chip_comparisons as comps
                           ON seeds.track = comps.track1
                          AND seeds.primary_id = comps.Query_ID
                      WHERE track = '%(track)s'
                      GROUP BY primary_id'''

        results = self.getDataFrame(statement)
      
        
        results["Similar motif found in"] = results["Similar motif found in"].apply(line2lines)
  
        img_tmp = '''.. image:: memechip.dir/%s_%i.png
   :scale: 25%%'''
        
        results["logo"] = results["ID"].apply(lambda x: img_tmp % (track,x))
        results["link"] = "`MEME-CHIP report <%s/index.html>`_" % resultsdir

        return results


class SequenceSubset(IntervalTracker):

    pattern = "(.+)_summary"

    def __call__(self, track):
        
        def _P(key):
            return P["%s_%s" % (track,key)]

        if _P("num_sequences"):
            which = str(_P("num_sequences"))
        elif _P("proportion"):
            if _P("min_sequences"):
                mins = _P("min_sequences")
            else:
                mins = 0
            which = "%s%% (minimum %i) of" % (str(_P("proportion") * 100),
                                           mins)
        else:
            which = "all"

        if _P("score") == "random":
            seqs = "%s sequences at random were selected" % which
        elif _P("score"):
            seqs = "The top %s sequences ranked by %s were selected" % (which, _P("score"))
        else:
            seqs = "%s sequences were selected" % which
        if _P("max_size"):
            seqs += ", up to a total of %ibp" % _P("max_size")

        if _P("halfwidth"):
            seqs += ".\n %ibp either side of the peak was used." % _P("halfwidth")
        else:
            seqs += ".\n The whole interval was used."

        if _P("masker"):
            seqs += "\n Sequences were masked with %s" % _P("masker")

        return seqs


class GeneSetComparision(IntervalTracker):

    def __call__(self, track, slice):

        track, prog = re.match("(.+)_([^_]+)", track).groups()
        table = self.slice2table[slice] % locals()
        column = self.slice2column[slice]

        return self.getValues("SELECT %(column)s FROM %(table)s")


class CompFullSubset(GeneSetComparision):

    slices = ["full", "subset"]
    pattern = "(.+)_motifseq_stats"
    slice2table = {"full": "%(track)s_composition",
                    "subset": "%(track)s_%(prog)s_motifseq_stats"}


class CompFullSubsetGC(CompFullSubset):

    slice2column = {"full": "CpG_density",
                    "subset": "pGC"}


class CompFullSubsetLength(CompFullSubset):

    slice2column = {"full": "end-start",
                    "subset": "(nA+nT+nG+nC +nN)"}


class CompPosNeg(GeneSetComparision):

    slices = ["positives", "negatives"]
    
    def getTracks(self):

        tracks = self.getValues('''SELECT track
                                   FROM %(prog)s_summary'''
                                )
        tracks = filter(lambda x: "_vs_" in x, tracks)
        return tracks

    def __call__(self, track, slice):

        groups = re.match("(.+)_vs_(.+)", track).groups()

        if slice == "positives":
            track = groups[0]
        else:
            track = groups[1]

        track = re.sub("-","_",track)
        return self.getValues(''' SELECT %(column)s
                                  FROM %(track)s_%(prog)s_motifseq_stats''')


class PosNegMemeGC(CompPosNeg):
    prog = "meme"
    column = "pGC"


class PosNegDremeGC(CompPosNeg):
    prog = "dreme"
    column = "pGC"


class PosNegMemeLength(CompPosNeg):
    prog = "meme"
    column = "nGC + nAT + nN"


class PosNegDremeLength(CompPosNeg):
    prog = "dreme"
    column = "nGC + nAT + nN"


class PosNegMemeN(CompPosNeg):
    prog = "meme"
    column = "pN"


class PosNegDremeN(CompPosNeg):
    prog = "dreme"
    column = "pN"


class TomTomResults(IntervalTracker):

    '''overview of tomtom results.'''

    pattern = "(.*)_tomtom$"

    def __call__(self, track, slice=None):

        data = self.getAll("""SELECT query_id, target_id, target_name,
        optimal_offset,pvalue,qvalue,evalue, overlap, query_consensus,
        target_consensus, orientation
        FROM %(track)s_tomtom""" % locals())

        resultsdir = os.path.abspath(
            os.path.join(EXPORTDIR, "tomtom", "%s.tomtom" % re.sub("_", "-", track)))
        if not os.path.exists(resultsdir):
            return []

        # format is: match_q_3_t_2_M01904
        # q_3: query_id
        # t_2: target database
        # M01904: target_id
        data['link'] = ["`tomtom <%s/tomtom.html#match_q_%s_t_2_%s>`_" % (resultsdir, target_id, target_name)
                        for target_id, target_name in zip(data['query_id'], data['target_id'])]
        return data


# class AnnotationsMatrix( DefaultTracker ):


#     def getSlices( self, subset = None ):
#         if subset: return subset
#         return []

#     def __call__(self, track, slice = None ):

#         result = odict()
#         rows = ("intergenic", "intronic", "upstream", "downstream", "utr", "cds", "other" )

#         statement = self.getStatement( slice )
#         data = self.get( statement % locals() )
#         levels = sorted(list(set( [ x[7] for x in data ] )))

#         for row in rows:
#             m = odict()
#             for l in levels: m[l] = 0
#             result[row] = m

#         map_level2col = dict( [(y,x) for x,y in enumerate(levels)] )
#         for intergenic, intronic, upstream, downstream, utr, coding, ambiguous, level in data:
#             col = level
#             for x,v in enumerate( (intergenic, intronic, upstream, downstream, utr, coding, ambiguous)):
#                 if v:
#                     row=rows[x]
#                     break
#             else:
#                 row = rows[-1]

#             result[row][col] += 1

#         return result

# class AnnotationsMotifs( AnnotationsMatrix ):
#     '''return a matrix with intervals stratified by motif presence
#     and location of the interval.
#     '''

#     mPattern = "_mast$"

#     def getStatement( self, slice = None ):

#         statement = '''
#         SELECT a.is_intergenic, a.is_intronic, a.is_upstream, a.is_downstream, a.is_utr, a.is_cds, a.is_ambiguous,
#           CASE WHEN m.nmatches > 0 THEN motif || '+' ELSE motif || '-' END
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id,
#         %(track)s_mast AS m ON m.id = i.interval_id'''

#         if slice is not None:
#             statement += " AND motif = '%(slice)s'"
#         return statement

# class AnnotationsPeakVal( AnnotationsMatrix ):
#     '''return a matrix with intervals stratified by peakval
#     and location of the interval.
#     '''
#     mPattern = "_annotations$"

#     def getStatement( self, slice = None ):

#         statement = '''
#         SELECT a.is_intergenic, a.is_intronic, a.is_upstream, a.is_downstream, a.is_utr, a.is_cds, a.is_ambiguous,
#         peakval
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id'''

#         return statement

# class AnnotationsPeakValData( DefaultTracker ):
#     '''return peakval for intervals falling into various regions.'''

#     def getSlices( self, subset = None ):
#         if subset: return subset
# return ("intergenic", "intronic", "upstream", "downstream", "utr",
# "cds", "other" )

#     def __call__(self, track, slice = None ):

#         if slice == "other": slice = "ambiguous"

#         statement = '''
#         SELECT peakval
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id AND a.is_%(slice)s ''' % locals()

#         return odict( (("peakval", self.getValues( statement )),) )
