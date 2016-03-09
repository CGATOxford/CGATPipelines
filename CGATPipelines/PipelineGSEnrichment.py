import pandas as pd
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
from CGATPipelines.Pipeline import cluster_runnable
import CGAT.IOTools as IOTools
import ast
import os
import rpy2
import rpy2.robjects as robjects
import rpy2.interactive as r
import rpy2.interactive.packages
r.packages.importr("hpar")
import scipy.stats as stats
import shutil

robjects.r('''
library("hpar")
hpaQuery <- function(tissue, level, supportive=T){
    data(hpaNormalTissue)
    intissue = grep(paste0("^", tissue, "*"),
                    ignore.case=T,
                    hpaNormalTissue$Tissue)

    if (level == "low"){str = "low|medium|high"}
    else if (level == "medium") {str = "medium|high"}
    else if (level == "high") {str = "high"}

    atlevel = grep(str, ignore.case=T,
                   hpaNormalTissue$Level)

    supp = grep("supportive", ignore.case=T,
                hpaNormalTissue$Reliability)

    all = intersect(intissue, atlevel)
    if (supportive == T) {
        all = intersect(all, supp)
    }
    return (as.vector(hpaNormalTissue[all,]$Gene))
}
''')


class AnnotationSet(object):
    '''
    Class for sets of annotations (not specific to the gene list).
    Generated for each set of annotations with a *_config.tsv ini
    file in the working dir.
    Stores mapping of genes to term IDs (currently only from a flat file)
    and the ontology context of each term (if a .obo file is provided).
    '''
    def __init__(self, optionsfile=None, annots_from_config=None, onts=dict(),
                 descs=dict()):
        self.optionsfile = optionsfile
        self.annots_from_config = annots_from_config
        self.onts = onts
        self.descs = descs
        self.prefix = str()
        self.GenesToTerms = dict()
        self.TermsToGenes = dict()
        self.TermsToOnt = dict()
        self.MappingFilePath = str()
        self.OntFilePath = str()
        self.Genes = set()

    def build(self):
        '''
        Reads and builds an annotation set from a mapping table and,
        if specified, an ontology file.
        '''
        if self.annots_from_config is True:
            self.AnnotationsFromConfig(self.optionsfile)
        elif self.annots_from_config is False:
            self.AnnotationsFromAnnotations(self.optionsfile)

        prefix = self.optionsfile.split("/")[-1]
        prefix = prefix.split("_")[0]
        self.prefix = prefix
        self.G2TFile = "%s_mapped.tsv" % prefix
        self.mapGenesAndTerms()
        self.Genes = set(self.GenesToTerms.keys())

        if "ontologies_%s" % prefix in self.onts:
            self.OntFilePath = self.onts["ontologies_%s" % prefix]
        else:
            self.OntFilePath = ""
        if "descs_%s" % prefix in self.descs:
            self.descfile = self.descs["termdescs_%s" % prefix]
        else:
            self.descfile = ""
        if os.path.exists(self.descfile):
            dDict = dict()
            id_desc = [x.strip().split("\t")
                       for x in open(self.descfile).readlines()]
            for pair in id_desc:
                dDict[pair[0]] = set([pair[1]])
            self.dDict = dDict
        else:
            self.dDict = dict()

        #  read ontology file if one is listed
        if self.OntFilePath is not None and len(self.OntFilePath) != 0:
            self.readOntFile()
            self.addChildren()
            self.reMap()

    def AnnotationsFromConfig(self, infile):
        opts = []
        for line in IOTools.openFile(self.optionsfile).readlines():
            line = line.replace(" ", "").strip()
            if len(line) != 0 and line[0] != "#":
                line = line.split("=")
                if len(line) != 1:
                    opts.append(line[1])
        if opts[2] != "0":
            statement = """
            python /ifs/devel/katherineb/cgat/scripts/annotations2annotations.py \
            -m standardise \
            -p %s \
            -l %s \
            --translate_ids %s \
            --translate_file %s \
            --translate_from %s \
            --translate_to %s \
            -k %s \
            -o %s \
            --delim1 "%s" \
            --delim2 "%s" \
            --ignore_delim_key %s \
            --ignore_delim_other %s \
            """ % tuple(opts)
        else:
            statement = """
            python /ifs/devel/katherineb/cgat/scripts/annotations2annotations.py \
            -m standardise \
            -p %s \
            -l %s \
            -k %s \
            -o %s \
            --delim1 "%s" \
            --delim2 "%s" \
            --ignore_delim_key %s \
            --ignore_delim_other %s \
            """ % tuple(opts[:2] + opts[6:])
        P.run()

    def AnnotationsFromAnnotations(self, mapped):
        GenesToTerms = self.unstowSetDict(mapped)
        TermsToGenes = dict()
        for gene in GenesToTerms:
            terms = GenesToTerms[gene]
            for term in terms:
                if term not in TermsToGenes:
                    TermsToGenes[term] = set()
                TermsToGenes[term].add(gene)
        self.GenesToTerms = GenesToTerms
        self.TermsToGenes = TermsToGenes
        shutil.copy(mapped, mapped.replace("mapping", "mapped"))

    def stow(self, outGenesToTerms, outTermsToGenes, outOnt=None,
             outGenes=None, noOnt=False):
        '''
        Writes the object to a set of output files so that it can be easily
        regenerated
        '''
        self.stowSetDict(self.GenesToTerms, outGenesToTerms)
        self.stowSetDict(self.TermsToGenes, outTermsToGenes)
        if outGenes is not None:
            self.stowGeneSet(outGenes)
        if self.OntFilePath is not None and noOnt is False:
            self.stowOnt(outOnt)
        elif noOnt is False:
            os.system("touch %s" % outOnt)
        else:
            pass

        dDict = self.dDict
        if len(dDict) != 0:
            self.stowSetDict(dDict, "%s_translate.tsv" % self.prefix)

    def rebuild(self, inGenesToTerms, inTermsToGenes, inOnt=None,
                inGenes=None):
        '''
        Regenerates the Annotation Set from the output files generated by
        stow().
        '''
        self.GenesToTerms = self.unstowSetDict(inGenesToTerms)
        self.TermsToGenes = self.unstowSetDict(inTermsToGenes)
        prefix = inGenesToTerms.split("/")[-1]
        prefix = prefix.split("_")[0]
        self.prefix = prefix
        if inOnt is not None:
            self.unstowOnt(inOnt)
        self.Genes = set(self.GenesToTerms.keys())
        if os.path.exists("%s_translate.tsv" % self.prefix):
            self.dDict = self.unstowSetDict("%s_translate.tsv" % self.prefix)
        else:
            self.dDict = dict()

    def mapGenesAndTerms(self):
        '''
        Builds a GenesToTerms dictionary based on the MappingDF from the
        mapping file.
        '''
        self.GenesToTerms = self.unstowSetDict(self.G2TFile)
        for key in self.GenesToTerms:
            terms = self.GenesToTerms[key]
            for term in terms:
                if term not in self.TermsToGenes:
                    self.TermsToGenes[term] = set()
                self.TermsToGenes[term].add(key)

    def reMap(self):
        TermsToOnt = self.TermsToOnt
        TermsToGenes = self.TermsToGenes
        Adict = dict()
        for term in TermsToGenes:
            Adict[term] = set()
            if term in TermsToOnt:
                A = TermsToOnt[term]
                A.getAll(TermsToOnt, "D")
                allT = set([A.ID]) | A.alldescendents
                for term2 in allT:
                    if term2 in TermsToGenes:
                        Adict[term] = Adict[term] | TermsToGenes[term2]

        Bdict = dict()
        for term in Adict:
            genes = Adict[term]
            for gene in genes:
                if gene not in Bdict:
                    Bdict[gene] = set()
                Bdict[gene].add(term)
        self.TermsToGenes = Adict
        self.GenesToTerms = Bdict

    def readOntFile(self):
        '''
        Reads an obo formatted file.
        Generates a dictionary - self.TermsToOnt
        '''
        TermsToOnt = self.TermsToOnt

        T = None
        s = 0
        with IOTools.openFile(self.OntFilePath) as infile:
            for line in infile:
                line = line.strip()
                if line.startswith("[Typedef]"):
                    s = 1
                if line.startswith("[Term]"):
                    T = ontologyTerm()
                if T is not None:
                    if s == 1:
                        if len(line) == 0:
                            s = 0
                    elif len(line) == 0:
                        TermsToOnt[T.ID] = T
                    else:
                        segs = line.split(": ")
                        if len(segs) >= 2:
                            T.add_var(segs[0], ":".join(segs[1:]))
        TermsToOnt[T.ID] = T
        self.TermsToOnt = TermsToOnt

    def addChildren(self):
        '''
        Adds "child" terms to each ID in a TermsToOnt attribute generated
        from an obo file.
        '''
        TermsToOnt = self.TermsToOnt
        for termid, term in self.TermsToOnt.items():
            parents = term.parents
            if parents is not None:
                for parent in parents:
                    self.TermsToOnt[parent].add_var("includes", termid)
        self.TermsToOnt = TermsToOnt

    def stowSetDict(self, adict, outfile):
        '''
        Stores a parsed mapping file in a flat file to load into the database
        and easily re-read later
        '''
        out = IOTools.openFile(outfile, "w")
        out.write("id1\tid2\n")

        for id1, id2 in adict.items():
            if len(id1) != 0:
                out.write("%s\t%s\n" % (id1,
                                        ",".join(id2)
                                        if len(id2) != 0 else "-"))
        out.close()

    def stowOnt(self, outfile):
        '''
        Stores a parsed ontology file in a flat file to load into the database
        and easily re-read later
        '''
        out = IOTools.openFile(outfile, "w")
        SB = self.TermsToOnt
        if len(SB.keys()) != 0:
            out.write("ID\tTerm\tDefinition\tParents\tChildren\n")
        for key in SB.keys():
            out.write("%s\t%s\t%s\t" % (SB[key].ID,
                                        SB[key].name,
                                        SB[key].defi))
            if len(SB[key].parents) == 0:
                out.write("-\t")
            else:
                out.write("%s\t" % (",".join(SB[key].parents)))
            if len(SB[key].children) == 0:
                out.write("-\n")
            else:
                out.write("%s\n" % (",".join(SB[key].children)))
        out.close()

    def unstowSetDict(self, infile):
        '''
        Regenerates a dictionary using a file generated by
        stowSetDict
        '''
        D = dict()
        i = 0
        with IOTools.openFile(infile) as inf:
            for line in inf:
                if i != 0:
                    line = line.strip().split("\t")
                    if line[0] not in D:
                        D[line[0]] = set()
                        for part in line[1].split(","):
                            D[line[0]].add(part)
                i += 1
        return D

    def unstowOnt(self, infile):
        '''
        Regenerates the TermsToOnt dictionary using the file generated by
        writeOntFlatFile.
        '''
        SB = self.TermsToOnt
        i = 0
        with IOTools.openFile(infile) as inf:
            for line in inf:
                if i != 0:
                    line = line.strip().split("\t")
                    T = ontologyTerm()
                    T.ID, T.name, T.defi = line[0:3]
                    T.parents = set(line[3].split(","))
                    T.children = set(line[4].split(","))
                    SB[T.ID] = T
                i += 1
        self.TermsToOnt = SB

    def subSetAnnotations(self, genelist):
        '''
        Removes all but a specified set of genes from the AnnotationSet
        '''
        GenesToTerms = self.GenesToTerms
        genes_subset = set(GenesToTerms.keys()) & genelist

        subGenesToTerms = dict()
        subTermsToGenes = dict()

        for gene in genes_subset:
            terms = GenesToTerms[gene]
            subGenesToTerms[gene] = terms
            for term in terms:
                if term not in subTermsToGenes:
                    subTermsToGenes[term] = set()
                subTermsToGenes[term].add(gene)

        self.TermsToGenes = subTermsToGenes
        self.GenesToTerms = subGenesToTerms

    def getGenesAnnotatedToTerm(self, term):
        if term in self.TermsToGenes:
            return self.TermsToGenes[term]
        else:
            return set()

    def getGenesNotAnnotatedToTerm(self, term):
        if term in self.TermsToGenes:
            return set(self.Genes - self.TermsToGenes[term])
        else:
            return set()


class ontologyTerm(object):
    '''
    Class for a term from an ontology file and associated information.
    Corresponds to one stanza from an obo file.
    '''
    def __init__(self):
        self.ID = None
        self.parents = set()
        self.name = None
        self.defi = None
        self.children = set()

    def add_var(self, nam, var):
        '''
        adds attributes to the ontologyTerm based on the id, is_a, name,
        def and includes tags in the obo file.
        '''
        if nam == "id":
            self.ID = var
        if nam == "is_a":
            var = var.split(" ")[0]
            self.parents.add(var)
        elif nam == "name":
            self.name = var
        elif nam == "def":
            self.defi = var
        elif nam == "includes":
            self.children.add(var)

    def getAll(self, SB, ad):
        '''
        Adds all the ancestors of a term (if ad = "A") or all the
        descendents of a term (if ad = "D") to the term as
        self.allancestors or self.alldescendents.  Takes time so only
        run when required.
        '''
        children = self.children
        parents = self.parents

        if ad == "A":
            ids = parents
        elif ad == "D":
            ids = children
        allids = set()
        done = set()

        while len(ids) != 0:
            allids = allids | ids
            for item in allids:
                if ad == "A":
                    newids = SB[item].parents
                elif ad == "D":
                    newids = SB[item].children
                ids = ids | newids
                done.add(item)
            ids = ids - done
        if ad == "A":
            self.allancestors = allids
        elif ad == "D":
            self.alldescendents = allids


class StatsTest(object):
    def __init__(self, foregroundfiles, backgroundfiles):
        foreground = AnnotationSet()
        background = AnnotationSet()
        if len(foregroundfiles) == 2:
            foregroundfiles.append(None)
        if len(backgroundfiles) == 2:
            backgroundfiles.append(None)
        foreground.rebuild(foregroundfiles[0], foregroundfiles[1],
                           foregroundfiles[2])
        background.rebuild(backgroundfiles[0], backgroundfiles[1],
                           backgroundfiles[2])
        self.foreground = foreground
        self.background = background

    def FisherExactTest(self, term):
        foreground = self.foreground
        background = self.background
        A = len(foreground.getGenesAnnotatedToTerm(term))
        B = len(foreground.getGenesNotAnnotatedToTerm(term))
        C = len(background.getGenesAnnotatedToTerm(term))
        D = len(background.getGenesNotAnnotatedToTerm(term))
        fishlist = ((A, B), (C, D))
        dDict = background.dDict
        if os.path.exists("%s_translate.tsv" % background.prefix):
            termtrans = str(dDict[term])
        else:
            termtrans = term
        return (stats.fisher_exact(fishlist), fishlist, termtrans)

    def TermForTerm(self):
        '''
        standard test for overrepresented GO term in foreground vs background
        FC = fold change
        '''
        allterms = self.background.TermsToGenes.keys()
        D = dict()
        for term in allterms:
            res = self.FisherExactTest(term)
            D[term] = res
        self.termscores = D

    def writeScores(self, outfile, bgnam):
        out = IOTools.openFile(outfile, "w")
        termscores = self.termscores
        n = len(termscores.keys())
        for pair in termscores.items():
            name = pair[0]
            result = pair[1]
            trans = result[2]
            pvalue = result[0][1]
            oddsratio = result[0][0]
            c1, c2 = str(result[1][0][0]), str(result[1][0][1])
            c3, c4 = str(result[1][1][0]), str(result[1][1][1])
            counts = "\t".join([c1, c2, c3, c4])
            bf_pvalue = pvalue * n
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                      % (bgnam, name, trans, counts, oddsratio,
                         pvalue, bf_pvalue))
        out.close()


def HPABackground(tissue, level, supportive=1):
    '''
    Generates a background set using human protein atlas annotations using the
    R hpar package.
    '''
    if int(supportive) == "1":
        supp = "T"
    else:
        supp = "F"
    genes = set(robjects.r.hpaQuery(tissue, level, supp))
    return genes


def translateGeneSet(genes, translatetab, fromcol, tocol):
    '''
    Translates a geneset from one set of identifiers to another using data
    in translatetab.
    fromcol = original ID type (column heading in translatetab)
    tocol = final ID (column heading in translatetab)
    '''
    i = 0
    tID = dict()
    with IOTools.openFile(translatetab) as inf:
        for line in inf:
            line = line.strip().split("\t")
            if i == 0:
                ind1 = line.index(fromcol)
                ind2 = line.index(tocol)

            else:
                if (len(line) > ind1) and (len(line[ind1]) != 0) and (
                        line[ind1] in genes):
                    if line[ind1] not in tID:
                        tID[line[ind1]] = set()
                    if len(line) > ind2:
                        tID[line[ind1]].add(line[ind2])

            i += 1
    tgenes = set()
    for gene in genes:
        if gene in tID:
            tgene = tID[gene]
            tgenes = tgenes | tgene
    return tgenes


def outputSet(alist, outfile):
    '''
    Writes a list of genes to an output file, one gene per line.
    '''

    out = IOTools.openFile(outfile, "w")
    for item in alist:
        out.write("%s\n" % item)
    out.close()


def inputSet(infile):
    '''
    Reads a set of genes from a file with one gene per line.
    '''
    return set(line.strip() for line in IOTools.openFile(infile).readlines())
