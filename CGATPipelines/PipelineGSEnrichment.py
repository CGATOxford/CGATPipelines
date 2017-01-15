import pandas as pd
import CGAT.IOTools as IOTools
import sqlite3
import os
import rpy2
import copy
import rpy2.robjects as robjects
import rpy2.interactive as r
import rpy2.interactive.packages
import scipy.stats as stats
from CGATPipelines.Pipeline import cluster_runnable
import CGAT.Experiment as E
import ast as ast
from toposort import toposort_flatten
r.packages.importr("hpar")


robjects.r('''
library("hpar")
# Uses the r hpar library to find human protein atlas terms
# meeting the criteria specified in the pipeline.ini

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


def removeNonAscii(s):
    '''
    Removes non-ascii characters from database terms (as some downloaded
    information has special characters which cause errors)
    '''
    return "".join(i for i in s if ord(i) < 128)


def getUnmapped(params):
    '''
    Parses the PARAMS['annotation_flatfiles'] strings from the pipeline.ini
    to find out which output files to create
    '''
    unmapped = dict()
    for p in params:
        if "_".join(p.split("_")[0:2]) == "annotation_flatfiles":
            L = [x.strip() for x in params[p].split("-")]
            for l in L:
                l = l.split(" ")
                if l[0] == "p":
                    pref = l[1]
                    unmapped[pref] = params[p]
    return unmapped


def getTables(dbname):
    '''
    Retrieves the names of all tables in the database.
    Groups tables into dictionaries by annotation
    '''
    dbh = sqlite3.connect(dbname)
    c = dbh.cursor()
    statement = "SELECT name FROM sqlite_master WHERE type='table'"
    c.execute(statement)
    tables = c.fetchall()
    c.close()
    dbh.close()
    D = {}
    for t in tables:
        tname = t[0].replace("ensemblg2", "").split("$")
        E.info(tname)
        ttype = tname[0]
        D.setdefault(ttype, [])
        D[ttype].append(tname[1])
    return D


def readDBTable(dbname, tablename):
    '''
    Reads the specified table from the specified database.
    Returns a list of tuples representing each row
    '''
    dbh = sqlite3.connect(dbname)
    c = dbh.cursor()
    statement = "SELECT * FROM %s" % tablename
    c.execute(statement)
    allresults = c.fetchall()
    c.close()
    dbh.close()
    return allresults


def getDBColumnNames(dbname, tablename):
    dbh = sqlite3.connect(dbname)
    res = pd.read_sql('SELECT * FROM %s' % tablename, dbh)
    dbh.close()
    return res.columns


@cluster_runnable
def translateGenelist(dbname, genelistfile, idtype):
    '''
    Translates a list of gene names from idtype to ensemblg based
    on the ensemblg2idtype$geneid database table.  This table needs
    to exist in the database.
    '''

    genelist = [line.strip() for line in
                IOTools.openFile(genelistfile).readlines()]
    trans = pd.DataFrame(
        readDBTable(dbname, "ensemblg2%s$geneid" % idtype))
    trans.columns = getDBColumnNames(dbname, "ensemblg2%s$geneid" % idtype)
    mergeon = set(trans.columns)
    mergeon.remove('ensemblg')
    mergeon = list(mergeon)[0]
    trans = trans[trans[mergeon].isin(genelist)]
    newgenelist = trans['ensemblg']
    return set(newgenelist.values)


@cluster_runnable
def untranslateGenelist(dbname, tab, genelistfile, idtype, outfile):
    genelist = [line.strip() for line in
                IOTools.openFile(genelistfile).readlines()]
    trans = pd.DataFrame(readDBTable(dbname, tab), columns=[idtype,
                                                            'ensemblg'])
    T = trans[idtype][trans['ensemblg'].isin(genelist)]
    T.to_csv(outfile, sep="\t", index=None)


@cluster_runnable
def writeList(genelist, outfile):
    '''
    Writes a list or set of genes to an output file, one per line
    '''
    outf = IOTools.openFile(outfile, "w")
    for gene in genelist:
        outf.write("%s\n" % gene)
    outf.close()


def getAllAncestorsDescendants(term, ontologyDict):
    '''
    Returns all the ancestors or descendents of a term at all levels of
    an ontology.
    ontologyDict should be a dictionary where keys are terms and values
    are sets of all the immediate children (for descendents) or parents
    (for ancestors) of those terms.
    '''
    ids = ontologyDict[term]
    allids = set()
    done = set()

    while len(ids) != 0:
        allids = allids | ids
        for item in allids:
            if item in ontologyDict:
                newids = ontologyDict[item]
            else:
                newids = set()
            ids = ids | newids
            done.add(item)
        ids = ids - done

    return allids


class AnnotationSet(object):
    '''
    Used to contain all the annotations associated with a particular
    annotation source.
    Consists of four dictionaries:
    1. GenesToTerms - keys are gene names, values are sets of terms mapped to
    these genes.

    2. TermsToGenes - keys are term names, values are sets of genes mapped to
    these terms (and their descendents if an ontology is provided)

    3. TermsToOnt - if an ontology file is provided, keys are terms, values are
    sets of terms which are direct parents of these terms (the is_a property
    in an obo or owl file)

    4. TermsToDetails -the annotation source usually provided other information
    besides a list of terms and the genes they are mapped to, e.g. each
    term will have a description.  TermsToDetails has a row for each
    term, the columns are any metadata about the term.

    An AnnotationSet also has a prefix used to store and regenerate the set,
    this is the stem of the output file names.
    e.g. prefix annotations.dir/hpo will generate
    annotations.dir/hpo_genestoterms.tsv
    annotations.dir/hpo_termstogenes.tsv
    annotations.dir/hpo_termstodetails.tsv
    annotations.dir/hpo_termstoont.tsv

    DetailsColumns contains the column names of the metadata in the
    TermsToDetails dictionary.  This is stored as the first row of the
    TermsToDetails output file.

    Methods allow the AnnotationSet to be:
    - Saved into a standard format (stow, stowSetDict, stowDetails)
    - Regenerated from this format (unstow, unstowSetDict, unstowDetails)
    - Reformatted (translateIDs, ontologise)
    '''

    def __init__(self, prefix):
        self.prefix = prefix
        self.GenesToTerms = dict()
        self.TermsToGenes = dict()
        self.TermsToOnt = dict()
        self.TermsToDetails = dict()
        self.DetailsColumns = []

    def unstow(self):
        '''
        Regenerates an AnnotationSet which has been stored as four
        files by stow().
        '''
        prefix = self.prefix
        self.GenesToTerms = self.unstowSetDict("%s_genestoterms.tsv" % prefix)
        self.TermsToGenes = self.unstowSetDict("%s_termstogenes.tsv" % prefix)
        self.TermsToOnt = self.unstowSetDict("%s_termstoont.tsv" % prefix)
        self.TermsToDetails, self.DetailsColumns = self.unstowDetails(
            "%s_termstodetails.tsv" % prefix)

    def stow(self, outprefix):
        '''
        Writes the object to a set of output files so that it can be easily
        regenerated
        '''
        outGenesToTerms = "%s_genestoterms.tsv" % outprefix
        outTermsToGenes = "%s_termstogenes.tsv" % outprefix
        outTermsToOnt = "%s_termstoont.tsv" % outprefix
        outTermsToDetails = "%s_termstodetails.tsv" % outprefix

        self.stowSetDict(self.GenesToTerms, outGenesToTerms, ["gene", "term"])
        self.stowSetDict(self.TermsToGenes, outTermsToGenes, ["term", "gene"])

        if self.TermsToOnt is not None:
            self.stowSetDict(self.TermsToOnt, outTermsToOnt, ['term', 'is_a'])
        else:
            os.system("touch %s" % outTermsToOnt)

        if self.TermsToDetails is not None:
            self.stowDetails(self.TermsToDetails, outTermsToDetails,
                             self.DetailsColumns)
        else:
            os.system("touch %s" % outTermsToDetails)

    def stowSetDict(self, adict, outfile, cnames):
        '''
        Stores a dictionary where values are sets or lists in a flat file with
        one column as the dictionary keys and the other a comma delimited
        list of the set of values for that key.
        e.g. {A:set("a", "b", "c")} would be stored as
        A    a,b,c
        '''
        out = IOTools.openFile(outfile, "w")
        out.write("%s\n" % ("\t".join(cnames)))

        for id1, id2 in list(adict.items()):
            out.write("%s\t%s\n" % (removeNonAscii(id1),
                                    removeNonAscii(",".join(id2))))
        out.close()

    def stowDetails(self, adict, outfile, cnames):
        '''
        Stores the TermsToDetails dictionary in a flat file
        This is slightly different to the other AnnotationSet dictionaries
        because order is important in the values in the dictionary so they
        are tuples rather than sets.
        The first line of the output file is the DetailsColumns column names
        '''
        out = IOTools.openFile(outfile, "w")
        out.write("%s\n" % ("\t".join(cnames)))
        for nam, val in list(adict.items()):
            tval = removeNonAscii("\t".join(val))
            out.write("%s\t%s\n" % (removeNonAscii(nam),
                                    tval))
        out.close()

    def unstowSetDict(self, infile):
        '''
        Regenerates a dictionary using a file generated by
        stowSetDict
        Reads a flat file where column one is dictionary keys and column two
        a comma delimited list of values for that key and generates a
        dictionary of sets.
        e.g. A    a,b,c would become  {A:set("a", "b", "c")}
        '''
        D = dict()
        i = 0
        with IOTools.openFile(infile) as inf:
            for line in inf:
                if i != 0:
                    line = line.strip().split("\t")
                    if line[0] not in D and len(line) > 1:
                        D[line[0]] = set()
                        for part in line[1].split(","):
                            D[line[0]].add(part)
                i += 1
        return D

    def unstowDetails(self, infile):
        '''
        Regenerates the dictionary stored by stowDetails.
        '''
        D = dict()
        i = 0
        with IOTools.openFile(infile) as inf:
            for line in inf:
                line = line.strip().split("\t")
                if i == 0:
                    cnames = line
                else:
                    D[line[0]] = tuple(line[1:])
                i += 1
        return D, cnames

    def translateIDs(self, dbname, idtype):
        '''
        Converts IDs throughout an AnnotationSet from idtype
        to ensemblg using information stored in the dbname database
        table ensemblg2idtype$geneid
        '''
        tab = readDBTable(dbname, "ensemblg2%s$geneid" % idtype)
        cols = getDBColumnNames(dbname, "ensemblg2%s$geneid" % idtype)
        mergeon = set(cols)
        mergeon.remove('ensemblg')
        mergeon = list(mergeon)[0]
        tabpd = pd.DataFrame(tab, columns=cols)

        genes = self.GenesToTerms
        tempdict = dict()
        for gene in self.GenesToTerms:
            tempdict[gene] = ",".join(genes[gene])
        df = pd.DataFrame.from_dict(tempdict, 'index')
        T = tabpd.merge(df, left_on=mergeon, right_index=True)
        T = T.drop(mergeon, 1)
        T.columns = ['ensemblg', 'terms']
        g2t = dict(list(zip(T['ensemblg'], T['terms'])))
        for key in g2t:
            g2t[key] = g2t[key].split(",")
        #  Replace the GenesToTerms and TermsToGenes dictionaries
        #  with copies with the translated gene names.  TermsToDetails
        #  and TermsToOnt don't contain gene names so they can stay the same.
        self.GenesToTerms = g2t
        self.TermsToGenes = self.reverseDict(self.GenesToTerms)

    def reverseDict(self, D):
        '''
        Inverts a dictionary of sets so that values are keys and
        keys are values
        e.g.
        {A: set(1, 2, 3) B: set(2, 3, 4)}
        -->
        {1: set(A), 2: set(A, B), 3: set(A, B), 4: set(B)}
        '''
        rD = dict()
        for key in D:
            res = D[key]
            for s in res:
                rD.setdefault(s, set())
                rD[s].add(key)
        return rD

    def ontologise(self):
        '''
        Takes the TermsToGenes and GenesToTerms dictionaries and
        corrects them for an AnnotationSet with a hierarchical ontology.
        If a gene is associated with a term, e.g. gene ENSG00000144061 is
        associated with the HPO term Nephropathy, it is also associated with
        all the ancestors of that term, e.g. ENSG00000144061 must also be
        associated with Abnormality of the Kidney.
        This function deals with this by taking the TermsToGenes
        dictionary and, for each term, taking the descendent terms,
        looking up their associated genes and adding them to the TermsToGenes
        set for the original term.

        '''
        TermsToGenes = self.TermsToGenes
        TermsToOntP = copy.copy(self.TermsToOnt)
        TermsToOntC = self.reverseDict(TermsToOntP)

        Adict = dict()
        # topologically sorts the terms in the ontology so that
        # every term is earlier in the list than all of its ancestors.
        sortedterms = toposort_flatten(self.TermsToOnt)
        sortedterms_p = []
        for s in sortedterms:
            if s in TermsToGenes:
                sortedterms_p.append(s)
        for term in sortedterms_p:
            Adict[term] = set()
            if term in TermsToOntC:
                # return descendents
                allids = getAllAncestorsDescendants(term, TermsToOntC)

                allids.add(term)
                for term2 in allids:
                    if term2 in TermsToGenes:
                        Adict[term] = Adict[term] | TermsToGenes[term2]

        self.TermsToGenes = Adict
        self.GenesToTerms = self.reverseDict(Adict)


class AnnotationParser(object):
    '''
    Base class for parsing an existing set of annotations into an
    AnnotationSet object
    '''

    def __init__(self, infile, prefix, options):
        self.infile = infile
        self.options = options
        self.AS = AnnotationSet(prefix)


class DBTableParser(AnnotationParser):
    '''
    Parses a set of database tables into an AnnotationSet

    3 table names are needed in the "options" dictionary:

    annot - each gene and each term it is annotated to.  There should be
    one gene and term per row, genes mapped to multiple terms are
    repeated multiple times e.g.
    gene term
    A    1
    A    2
    B    2
    details - each term and its metadata

    ont - each term and its parent terms (in an ontology) - can be None if the
    annotation is not hierarchical.

    prefix is the prefix for the output files, so files will be named
    prefix_genestoterms.tsv
    prefix_termstogenes.tsv
    prefix_termstoont.tsv
    prefix_termstodetails.tsv
    '''

    def __init__(self, infile, prefix, options):
        AnnotationParser.__init__(self, infile, prefix, options)
        self.annot = options['annot']
        self.details = options['details']
        self.ont = options['ont']

    def run(self):
        '''
        Runs the functions to convert the database into an AnnotationSet.
        '''
        self.AS.GenesToTerms, self.AS.TermsToGenes = self.annotsToDicts()
        if self.details is not None:
            (self.AS.TermsToDetails,
             self.AS.DetailsColumns) = self.detailsToDict()
        else:
            self.AS.TermsToDetails = None
        if self.ont is not None:
            self.AS.TermsToOnt = self.ontToDict()
            self.AS.ontologise()

        else:
            self.AS.TermsToOnt = None

    def annotsToDicts(self):
        '''
        Takes the annot table and builds two dictionaries to represent it.
        In rdict1 each gene is a key and each value is a set listing the
        terms associated with that gene.
        In rdict2 each term is a key and each value is a set listing the
        genes associated with that term.
        These become the GenesToTerms and TermsToGenes dictionaries of the
        AnnotationSet.
        '''
        allresults = readDBTable(self.infile, self.annot)
        cnames = getDBColumnNames(self.infile, self.annot)

        if cnames[0] == "ensemblg":
            ind0 = 0
            ind1 = 1
        else:
            ind0 = 1
            ind1 = 0
        rdict1 = {}
        rdict2 = {}

        for result in allresults:
            rdict1.setdefault(result[ind0], set())
            rdict1[result[ind0]].add(result[ind1])

            rdict2.setdefault(result[ind1], set())
            rdict2[result[ind1]].add(result[ind0])

        return rdict1, rdict2

    def detailsToDict(self):
        '''
        Takes the details table and builds a dictionary to represent it.
        Each dictionary value is a list corresponding to a table row
        (this needs to be ordered) and contains "-" where a column is blank.
        Each dictionary key is a term ID.
        '''
        allresults = readDBTable(self.infile, self.details)
        D = {}
        for line in allresults:
            res = []
            for x in line[1:]:
                if x is None:
                    res.append("-")
                else:
                    res.append(x)
            D[line[0]] = res

        cnames = getDBColumnNames(self.infile, self.details)
        return D, cnames

    def ontToDict(self):
        '''
        Reads the ontology table from the database into a dictionary
        Each dictionary key is a term and each value a set containing the
        direct parents of that term (the is_a attribute of the ontology)
        '''
        allresults = readDBTable(self.infile, self.ont)
        D = {}
        for result in allresults:
            D.setdefault(result[0], set())
            D[result[0]].add(result[1])
        return D


class FlatFileParser(AnnotationParser):
    '''
    Used to add annotations from flat files outside of the database
    and parse them with the annotations from the database.
    Allows the user to add custom annotations.
    '''

    def __init__(self, instring):
        '''
        Reads options strings as specfied in self.readOptions()
        '''
        options = self.readOptions(instring)
        AnnotationParser.__init__(self, options['l'],
                                  options['p'], options)

    def run(self, dbname):
        '''
        Reads the flat file into the four AnnotationParser dictionaries.
        '''
        self.readFile()
        self.AS.TermsToGenes = self.AS.reverseDict(self.AS.GenesToTerms)
        if self.options['ont'] is not None:
            self.AS.TermsToOnt = self.readOntFile()
            self.AS.ontologise()
        if self.options['ids'] != "ensemblg":
            self.AS.translateIDs(dbname, self.options['ids'])

    def readOptions(self, optionsstring):
        '''
        Reads an "optionsstring" from the pipeline.ini file and parses
        this into a dictionary.
        The options string needs to contain these options:
        -p -  prefix for output files

        -l - path to the flat file

        -k - this can be either the column number (0 indexed) or the column
        heading for the column containing the names of the genes

        -o - this can be either the column number or the column heading for
        the column containing the term IDs

        These options are optional:
        -d1 - this is the column delimiter
        If nothing is specified "\t" or tab delimited will be assumed.
        Otherwise give the option in double quotes e.g. -d1 "|"

        -d2 - the within column delimiter used when multiple genes and terms
        are specified in the same table row
        e.g. gene1,gene2,gene3   term1,term2 would need
        a -d2 parameter of ","

        -idk - if this is 1, ignore the column delimiter in the column
        containing the gene names
        -ido - if this is 1, ignore column delimiter in the column
        sometimes a within column delimiter is used not as a delimiter in
        another column e.g. in the line A1,A2,A3    Kidney, Abnormality
        of - you may want to
        seperate A1,A2,A3 based on the comma but ignore the comma in
        Kidney, Abnormality of.  These parameters allow you to ignore the
        delimiter in the gene or the term column. containing the term names

        -ont - path to ontology hierarchy
        -ids - type of ID, e.g. symbol, entrez, ensemblg - default is ensemblg
        '''
        options = dict()
        sections = optionsstring.strip().split("-")
        for section in sections:
            section = section.strip()
            if len(section) != 0:
                t, v = section.split(" ")
                options[t] = v
        options.setdefault('l', None)
        options.setdefault('k', '0')
        options.setdefault('o', '1')
        options.setdefault('d1', "\t")
        options.setdefault('d2', None)
        options.setdefault('idk', 0)
        options.setdefault('ido', 0)
        options.setdefault('ont', None)
        options.setdefault('ids', 'ensemblg')

        options['d1'] = ast.literal_eval("'%s'" % options['d1'])
        options['d2'] = ast.literal_eval("'%s'" % options['d2'])

        return options

    def readFile(self):
        '''
        Reads the input files to memory.
        '''
        # Determines parameters depending if columns are specified with
        # integers or with column names
        isd = False
        options = self.options
        if options['k'].isdigit():
            usecols = [int(s) for s in [options['k'], options['o']]]
            isd = True

        else:
            allcols = IOTools.openFile(
                options['l']).readline().split(options['d1'])
            k = allcols.index(options['k'])
            o = allcols.index(options['o'])
            usecols = [k, o]

        inf = IOTools.openFile(options['l']).readlines()
        (self.AS.GenesToTerms, self.AS.TermsToDetails,
         self.AS.DetailsColumns) = self.makeDict(inf, usecols, isd)

    def makeDict(self, inf, usecols, isd):
        '''
        Builds dictionaries based on the input files
        Dictionary keys are strings and values are sets.
        '''
        delim1 = self.options['d1']
        delim2 = self.options['d2']
        ignore1 = self.options['idk']
        ignore2 = self.options['ido']

        D = dict()
        detaildict = dict()
        i = 0
        for line in inf:
            line = line.strip().split(delim1)
            idpart = line[usecols[0]]
            termpart = line[usecols[1]]
            # if at line 0 and there are column headings
            if i == 0 and isd is False:
                cnames = copy.copy(line)
                cnames.remove(idpart)
                i += 1
                continue

            if delim2 is not None:
                if ignore1 == 0:
                    idparts = idpart.split(delim2)
                else:
                    idparts = [idpart]
                if ignore2 == 0:
                    termparts = termpart.split(delim2)
                else:
                    termparts = [termpart]
            else:
                idparts = [idpart]
                termparts = [termpart]
            # metadata needs to be stored in order
            for id in idparts:
                id = id.replace(" ", "")
                D.setdefault(id, set())
                for term in termparts:
                    term = term.replace(" ", "_")
                    term = term.replace(",", "")
                    D[id].add(term)
                    if term not in detaildict:
                        w = copy.copy(line)
                        w.remove(idpart)
                        w.remove(termpart)
                        w = tuple(w)
                        detaildict[term] = w
        # add default column names if none in the file
        if isd is True:
            cnames = tuple(["c_%s" % j for j in range(len(line) - 1)])
        return D, detaildict, cnames

    def readOntFile(self):
        '''
        Reads an obo formatted file.
        Generates a dictionary - self.TermsToOnt
        Keys are terms, values are all immediate parents of that term - the
        "is_a" argument in the obo file.
        '''
        TermsToOnt = dict()
        Tname = None
        isas = set()
        with IOTools.openFile(self.options['ont']) as infile:
            for line in infile:
                line = line.strip()
                if line.startswith("id"):
                    if Tname is not None:
                        TermsToOnt[Tname] = isas
                        isas = set()
                    Tname = line.split(": ")[1]
                elif line.startswith("is_a"):
                    isas.add(":".join(line.split(": ")[1:]).split(" ! ")[0])
        TermsToOnt[Tname] = isas
        return TermsToOnt


class EnrichmentTester(object):
    '''
    Runs statistical tests for enrichment in the "foreground" gene list
    compared to the "background" gene list for all terms in an AnnotationSet
    which are mapped to one or more genes in the foreground list.
    Outputs these results to a file in a standard format.
    '''

    def __init__(self, foreground, background, AS, runtype,
                 testtype, correction, thresh, outfile, outfile2,
                 idtype, dbname):
        allgenes = set(AS.GenesToTerms.keys())

        #  read the list of background genes, remove genes not in the
        #  AnnotationSet
        self.background = (set([line.strip()
                                for line in
                                IOTools.openFile(background).readlines()]) &
                           allgenes)

        # read the list of background genes, remove genes not in the
        # AnnotationSet and genes not in the background (all genes in
        # the foreground should also be in the background)
        self.foreground = (set([line.strip()
                                for line in
                                IOTools.openFile(foreground).readlines()]) &
                           allgenes) & self.background

        original_fg = set([line.strip() for line in
                           IOTools.openFile(foreground.replace(
                               "clean_", "")).readlines()])
        if os.path.exists(background.replace("clean_", "")):
            original_bg = set([line.strip() for line in
                               IOTools.openFile(foreground.replace(
                                   "clean_", "")).readlines()])
        else:
            original_bg = None
        self.ofg = original_fg
        self.obg = original_bg
        self.AS = AS
        self.runtype = runtype
        self.testtype = testtype
        self.correction = correction
        self.thresh = thresh
        terms = set()
        # Only terms which are associated with a gene in the foreground
        # are tested
        for gene in self.foreground:
            terms = terms | AS.GenesToTerms[gene]
        self.terms = terms
        self.outfile = outfile
        self.outfile2 = outfile2
        self.idtype = idtype
        self.dbname = dbname

    def writeStats(self, resultsdict, outfile, outfile2, writegenes,
                   host, ngenes):
        '''
        Writes the stats test results to an output table
        Each row has all the detail from the TermsToDetails table appended
        to it.  If the original table this was generated from did not have
        column names, the columns are named c1 to cn.
        outfile2 has the same information but only for significant tests
        and sorted by oddsratio.
        '''
        parsedresults = []
        i = 0
        tdict_bg = dict()
        tdict_fg = dict()
        for term in resultsdict:
            res = resultsdict[term]
            # odds ratio
            OR = round(res[0], 4)
            # p value
            p = round(res[1], 4)
            # adjusted p value
            padj = round(res[2], 4)
            # significance
            sig = "*" if res[3] is True else "-"
            # n genes associated with and not associated with
            # this term in the foreground
            A, B = res[4][0]
            # n genes associated with and not associated with
            # this term in the background
            C, D = res[4][1]

            L = [term, str(A), str(B), str(C), str(D),
                 str(OR), str(p), str(padj), sig]

            # add metadata
            L += self.AS.TermsToDetails[term]
            parsedresults.append(L)
            tdict_fg[str(term)] = res[-2]
            tdict_bg[str(term)] = res[-1]
            i += 1

        cols = ['term_id', 'fg_genes_mapped_to_term',
                'fg_genes_not_mapped_to_term',
                'bg_genes_mapped_to_term',
                'bg_genes_not_mapped_to_term', 'odds_ratio',
                'pvalue', 'padj', 'significant'] + self.AS.DetailsColumns[1:]
        df = pd.DataFrame(parsedresults, columns=cols)
        df = df.sort_values('padj')
        df.to_csv(outfile, sep="\t", index=False)

        df2 = df[df['significant'] == "*"]
        df2 = df2.sort_values('odds_ratio', ascending=False)
        df2.to_csv(outfile2, sep="\t", index=False)

        if writegenes == 1:
            #   writes the gene IDs annotated to the n most significantly
            #   enriched terms in this database to tsv files
            terms = df2['term_id'][0: ngenes]
            if len(terms) > 0:
                for term in terms:
                    fg = outfile.replace(".tsv", "_fg_genes_%s.tsv"
                                         % (term.replace(":", "_")))
                    fgo = IOTools.openFile(fg, "w")
                    bg = outfile.replace(".tsv", "_bg_genes_%s.tsv"
                                         % (term.replace(":", "_")))
                    bgo = IOTools.openFile(bg, "w")
                    fggenes = tdict_fg[term]
                    bggenes = tdict_bg[term]

                    for gene in fggenes:
                        fgo.write("%s\n" % gene)
                    for gene in bggenes:
                        bgo.write("%s\n" % gene)
                    fgo.close()
                    bgo.close()


class TermByTermET(EnrichmentTester):
    '''
    Standard method to look for enrichment - test each term
    individually.
    '''

    def __init__(self, foreground, background, AS, runtype,
                 testtype, correction, thresh, outfile, outfile2, idtype,
                 dbname):

        EnrichmentTester.__init__(self, foreground, background, AS, runtype,
                                  testtype, correction, thresh, outfile,
                                  outfile2, idtype, dbname)

    def run(self, writegenes, host, ngenes):
        results = dict()
        for term in self.terms:
            GenesWith = self.AS.TermsToGenes[term]
            GenesWithout = set(self.AS.GenesToTerms.keys()) - GenesWith

            if self.testtype == "Fisher":
                ST = FisherExactTest(term,
                                     self.foreground,
                                     self.background,
                                     GenesWith,
                                     GenesWithout,
                                     self.correction,
                                     self.thresh,
                                     len(self.terms), self.idtype,
                                     self.dbname,
                                     self.ofg,
                                     self.obg)
                res = ST.run()
                results[term] = res
        self.writeStats(results, self.outfile, self.outfile2,
                        writegenes, host, ngenes)


class EliminateET(EnrichmentTester):
    '''
    Implementation of the "elim" method of looking for enriched groups
    described in
    http://bioinformatics.oxfordjournals.org/content/22/13/1600.full.pdf
    Starts at the bottom of an ontology and works up towards the root - if
    a term at the bottom of the tree is significantly enriched
    in a geneset, its ancestors are eliminated from the analysis
    as these will also be enrihced but are less informative.
    '''

    def __init__(self, foreground, background, AS, runtype, testtype,
                 correction, thresh, outfile, outfile2, idtype, dbname):
        EnrichmentTester.__init__(self, foreground, background, AS, runtype,
                                  testtype, correction, thresh, outfile,
                                  outfile2, idtype, dbname)

    def run(self, writegenes, host, ngenes):
        TermsToOntP = copy.copy(self.AS.TermsToOnt)
        TermsToOntC = self.AS.reverseDict(TermsToOntP)

        # sort terms into "topological" order - every term comes earlier in the
        # list than all of its descendents
        sortedL = toposort_flatten(TermsToOntC)[::-1]

        # exclude terms not associated with any foreground genes
        p_sortedL = []
        for s in sortedL:
            if s in self.terms:
                p_sortedL.append(s)

        # if no terms associated with foreground function output is empty
        if not p_sortedL:
            results = dict()
            self.writeStats(results, self.outfile, self.outfile2,
                            writegenes, host, ngenes)
            return

        # find the "level" of each gene in the ontology - the number of
        # steps in the longest path from the node to the root
        levels = dict()
        results = dict()
        for term in p_sortedL:
            levels.setdefault(term, 0)
            if term in TermsToOntP:
                # what are the immediate parents of this term
                parents = TermsToOntP[term]
                maxi = 0
                if len(parents) != 0:
                    for parent in parents:
                        levels.setdefault(parent, 0)
                        if levels[parent] > maxi:
                            maxi = levels[parent]
                # the level of the term is 1 more than the maximum level of
                # its parents
                levels[term] = levels[term] + maxi + 1

        # invert the "levels" dictionary to give a dictionary where keys
        # are levels (integers) and values are sets of terms at that
        # level
        bylevel = dict()
        for term in levels:
            if term in self.terms:
                level = levels[term]
                bylevel.setdefault(level, set())
                bylevel[level].add(term)

        # find the highest level in the ontology
        maxLevel = max(bylevel.keys())
        markedGenes = dict()
        allgenes = set(self.AS.GenesToTerms.keys())
        for term in self.terms:
            markedGenes[term] = set()

        # iterate through the levels starting at the highest level - the
        # bottom of the tree
        for j in range(1, maxLevel)[::-1]:
            for term in bylevel[j]:
                genes = self.AS.TermsToGenes[term]

                # remove "marked genes" from the list of genes associated
                # with the term - these are genes which have already been
                # associated with a descendent of the term
                GenesWith = genes - markedGenes[term]
                GenesWithout = allgenes - genes
                if self.testtype == "Fisher":
                    ST = FisherExactTest(term, self.foreground,
                                         self.background,
                                         GenesWith,
                                         GenesWithout,
                                         self.correction,
                                         self.thresh,
                                         len(self.terms), self.idtype,
                                         self.dbname, self.ofg, self.obg)
                    R = ST.run()
                    results[term] = R

                if R[3] is True:
                    # "mark" the genes mapped to this term so
                    # they are not also mapped to ancestors of the term
                    ancs = getAllAncestorsDescendants(term, TermsToOntP)
                    for anc in ancs:
                        if anc in self.terms:
                            markedGenes[anc] = markedGenes[anc] | GenesWith
        self.writeStats(results, self.outfile, self.outfile2,
                        writegenes, host, ngenes)


class StatsTest(object):
    '''
    Container for StatsTest objects corresponding to different types
    of statistical test for enrichment.
    '''

    def __init__(self, term, foreground, background, GenesWith, GenesWithout,
                 correction, thresh, ntests, idtype, dbname, ofg, obg):
        self.term = term
        self.foreground = foreground
        self.background = background
        # Genes which are associated with the term
        self.GenesWith = GenesWith
        # Genes which are not associated with the term
        self.GenesWithout = GenesWithout
        self.correction = correction
        self.thresh = thresh
        self.ntests = ntests
        self.idtype = idtype
        self.dbname = dbname
        self.ofg = ofg
        self.obg = obg
        self.collapse()

    def collapse(self):
        # collapses the lists of genes back to the original id type
        # as different types of ID do not have a 1:1 relationship
        # avoids this if idtype already ensemblg
        if self.idtype != "ensemblg":
            db = sqlite3.connect(self.dbname)
            tab = pd.read_sql_query(
                "SELECT * FROM ensemblg2%s$geneid" % self.idtype, db)
            id = list(tab.columns)
            id.remove('ensemblg')
            id = id[0]
            self.GenesWith = set(tab[id][tab['ensemblg'].isin(self.GenesWith)])
            self.GenesWithout = set(
                tab[id][tab['ensemblg'].isin(self.GenesWithout)])
            self.foreground = set(
                tab[id][tab['ensemblg'].isin(self.foreground)])
            self.background = set(
                tab[id][tab['ensemblg'].isin(self.background)])

            self.foreground = self.foreground & self.ofg
            if self.obg is not None:
                self.background = self.background & self.obg

    def correct(self, pvalue):
        '''
        Correction for multiple testing
        '''
        # Bonferroni correction
        if self.correction == "bon":
            return pvalue * self.ntests

    def sig(self, pvalue):
        '''
        Test for significance
        '''
        if pvalue <= self.thresh:
            return True


class FisherExactTest(StatsTest):
    '''
    Runs a Fisher Exact Test (as implemented in scipy.stats)
    on a term.
    Uses a 2x2 contingency table with counts of:
    A - Genes in foreground annotated to term
    B - Genes in foreground not annotated to term
    C - Genes in background annotated to term
    D - Genes in background not annotated to term

    '''

    def __init___(self, term, foreground, background, GenesWith, GenesWithout,
                  correction, thresh, ntests, idtype, dbname):
        StatsTest.__init__(self, term, foreground, background, GenesWith,
                           GenesWithout, correction, thresh, ntests, idtype,
                           dbname)

    def run(self):
        A = self.GenesWith & self.foreground
        B = self.GenesWithout & self.foreground
        C = self.GenesWith & self.background
        D = self.GenesWithout & self.background

        fishlist = ((len(A), len(B)), (len(C), len(D)))

        FT = stats.fisher_exact(fishlist)
        OR, p = FT
        padj = self.correct(p)
        if padj > 1:
            padj = 1
        significant = self.sig(padj)

        return OR, p, padj, significant, fishlist, A, C

# functions below here correspond to specific steps in the
# pipeline_enrichment pipeline - they are written as functions
# so the cluster_runnable decorater can be used.


@cluster_runnable
def getDBAnnotations(infile, outfiles, dbname):
    '''
    Retrieves the database tables and runs DBTableParser as
    appropriate to generate AnnotationSets.
    '''
    dbtabs = getTables(dbname)
    for tab in dbtabs:
        if 'annot' in dbtabs[tab]:
            stem = tab
            Tdict = {'stem': stem,
                     'annot': "ensemblg2%s$annot" % stem,
                     'details': None,
                     'ont': None}
            if 'ont' in dbtabs[tab]:
                Tdict['ont'] = "%s$ont" % stem
            if 'details' in dbtabs[tab]:
                Tdict['details'] = "%s$details" % stem
            D = DBTableParser(dbname, stem, Tdict)
            D.run()
            D.AS.stow("annotations.dir/%s" % stem)


@cluster_runnable
def getFlatFileAnnotations(optionsstring, outstem, dbname):
    '''
    Retrieves annotations from a flat file.
    '''
    D = FlatFileParser(optionsstring)
    D.run(dbname)
    D.AS.stow(outstem)


@cluster_runnable
def translateAnnotations(infiles, outfiles, dbname):
    '''
    Translate a set of user specified annotations to use
    ensemblg identifiers, so all annotations have a consistent
    identifier type.
    '''
    parts = infiles[0].split("_")
    # regenerate the stored AnnotationSet
    AS = AnnotationSet("_".join(parts[0:3]))
    AS.unstow()
    AS.translateIDs(dbname, parts[2])
    AS.stow(parts[1])


@cluster_runnable
def cleanGeneLists(infile, outfile, idtype, dbname):
    '''
    Cleans a list of genes - removes duplicates and translates
    from idtype to ensemblg
    '''

    if idtype == "ensemblg":
        cleangenes = set([line.strip()
                          for line in IOTools.openFile(infile).readlines()])
    else:
        cleangenes = translateGenelist(dbname, infile, idtype)
    writeList(cleangenes, outfile)


@cluster_runnable
def HPABackground(tissue, level, supportive, outfile):
    '''
    Generates a background set using human protein atlas annotations using the
    R hpar package and thresholds set in pipeline.ini
    '''
    if int(supportive) == "1":
        supp = "T"
    else:
        supp = "F"
    genes = set(robjects.r.hpaQuery(tissue, level, supp))
    writeList(genes, outfile)


@cluster_runnable
def foregroundsVsBackgrounds(infiles, outfile, outfile2,
                             testtype, runtype, correction, thresh, dbname,
                             writegenes, host, ngenes, idtype):
    '''
    Runs the appropriate EnrichmentTester method according to
    the parameters specified in the pipeline.ini.
    '''
    annots = infiles[-1].replace("_genestoterms.tsv", "")

    # Regenerate the appropriate AnnotationSet
    AS = AnnotationSet(annots)
    AS.unstow()
    foreground = infiles[1]
    background = infiles[0]

    # Elim method doesn't work if there is no ontology - term by
    # term is used instead
    if runtype == "elim":
        if len(AS.TermsToOnt) == 0:
            runtype = "termbyterm"

    if runtype == "termbyterm":
        ET = TermByTermET(foreground, background, AS,
                          runtype, testtype, correction, thresh,
                          outfile, outfile2, idtype, dbname)
    elif runtype == "elim":
        ET = EliminateET(foreground, background, AS,
                         runtype, testtype, correction, thresh,
                         outfile, outfile2, idtype, dbname)
    ET.run(writegenes, host, ngenes)
