import CGATPipelines.Pipeline as P
from Bio import Entrez
import numpy as np
import httplib2
import json as json
import sqlite3
from intermine.webservice import Service
import string
import re
import os
import datetime
import xml.etree.ElementTree as ET
import pandas as pd
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import urllib2 as urllib2
from CGATPipelines.Pipeline import cluster_runnable

def removeNonAscii(s):
    '''
    Removes non-ascii characters from database terms (as some downloaded
    information has special characters which cause errors)
    '''
    if type(s) is str:
        return "".join(i for i in s if ord(i) < 128)
    else:
        return s


def readGeneList(filename):
    return [line.strip() for line in IOTools.openFile(filename).readlines()]


def getSymbols(hostfile):
    return [line.strip().split("\t")[0]
            for line in IOTools.openFile(hostfile).readlines()]


class APIAnnotation(object):
    def __init__(self, prefix, datasource, outdb, options=None, ohost=None,
                 email=None):
        self.prefix = prefix
        self.datasource = datasource
        self.outdb = outdb
        self.options = options
        self.ohost = ohost
        self.email = email
        self.resultsz = dict()
        self.resultspd = dict()

    def runall(self, genelist, fields=None, load=True):
        T = self.test()
        assert T == 1, "API test failed for %s" % self.prefix
        self.download(genelist, fields)
        self.parse()
        if load is True:
            self.loadPDTables()
            self.loadZippedTables()

    def translate(self, dataframe, tfrom, tto):
        assert ((os.path.exists('%s2%s.tsv' % (tfrom, tto))) |
                (os.path.exists('%s2%s.tsv' % (tto, tfrom)))), """
                translation from %s to %s needs to be available prior
                to %s parsing""" % (tfrom, tto, self.prefix)
        try:
            ensdf = self.getTranslation('%s2%s.tsv' % (tfrom, tto))
        except:
            ensdf = self.getTranslation('%s2%s.tsv' % (tto, tfrom))

        dataframe[tfrom] = dataframe[tfrom].astype(str)
        ensdf[tfrom] = ensdf[tfrom].astype(str)
        return dataframe.merge(
            ensdf, left_on=tfrom, right_on=tfrom).drop(tfrom, 1)

    def loadZippedTables(self, PS=False):
        '''
        Adds a database table, named tablename, to the database dbname.
        Cnames - column names as a Python list
        Ctypes - column types as strings in a Python list
        e.g. ['varchar (250)', 'varchar(25)', 'int']
        zipped - a list of tuples to fill the table, each item in the list
        should correspond to a row and should be a tuple of the value for
        each column in that row

        E.g.
        add_db_table(
        "csvdb", "entrez2symbol",
        ["entrez_id", "gene_symbol"],
        ['int', 'varchar (25)'],
        [(1111, 'aaa'), (2222, 'bbb'), (3333, 'ccc)])

        would create the table "entrez2symbol" in the database "csvdb" with the
        layout
        entrez_id    gene_symbol
        1111         aaa
        2222         bbb
        3333         ccc

        If PS is True SQL statements are printed before they are run - used for
        debugging.
        '''
        dbname = self.outdb
        # connect to the database
        dbh = sqlite3.connect(dbname)
        c = dbh.cursor()

        for tabkey in self.resultsz:
            tab = self.resultsz[tabkey]
            tablename = tabkey
            cnames = tab[1]
            ctypes = tab[2]
            zipped = tab[0]
            outf = (IOTools.openFile(tabkey + ".load", "w"))
            #  check there are the right number of column names for the dataset
            assert len(cnames) == len(
                zipped[0]), "length of column names doesn't \
                match number of fields"

            # build parts of the SQL statement
            cnamestring = ",".join(cnames)
            cnamestypes = ",".join([" ".join(z) for z in zip(cnames, ctypes)])
            qs = ",".join("?" * len(cnames))

            # create the table with the appropriate columns and column types
            statement = 'CREATE TABLE IF NOT EXISTS %s (%s);' % (tablename, cnamestypes)
            c.execute(statement)

            # populate the table with the data
            statement = '''INSERT INTO %s (%s) VALUES (%s);''' % (tablename,
                                                                  cnamestring,
                                                                  qs)
            outf.write("%s\n" % statement)
            if PS is True:
                E.info(statement)
            c.executemany(statement, zipped)
            outf.write("%s\n" % statement)
            dbh.commit()
            outf.close()
        # disconnect
        c.close()
        dbh.close()

    def loadPDTables(self):
        dbh = sqlite3.connect(self.outdb)
        for tabkey in self.resultspd:
            tab = self.resultspd[tabkey]
            tab.to_sql(tabkey, dbh, index=False, if_exists='replace')
            try:
                tab.to_csv("%s.load" % tabkey, sep="\t", index=False)
            except:
                for col in tab.columns:
                    for ind in tab.index:
                        tab.loc[ind, col] = removeNonAscii(tab.loc[ind, col])
                tab.to_csv("%s.load" % tabkey, sep="\t", index=False)
        dbh.close()

    def storeTranslation(self, genelist1, genelist2, cnames, out):
        df = pd.DataFrame(zip(genelist1, genelist2), columns=cnames)
        df.to_csv(out, sep="\t", index=False)

    def getTranslation(self, translation):
        df = pd.read_csv(translation, sep="\t")
        return df


class EntrezAnnotation(APIAnnotation):
    def __init__(self, prefix, outdb, email):
        APIAnnotation.__init__(self, prefix, "Entrez", outdb)
        Entrez.email = email

    def test(self):
        '''
        Test Entrez API - Look up symbol APOBEC3G, Entrez ID 60489 should be
        amongst the results
            '''
        E = Entrez.esearch(db="gene", term="(APOBEC3G[Preferred+Symbol])",
                           retmode="text", retmax=1000000)
        res = Entrez.read(E)
        E.close()
        if "60489" in res['IdList']:
            return 1
        else:
            return 0


class EntrezGeneAnnotation(EntrezAnnotation):
    def __init__(self, outdb, email):
        EntrezAnnotation.__init__(self, "entrezgene", outdb, email)

    def download_all(self, host):
        '''
        Gets all the Gene IDs for a particular host, specified in PARAMS, from
        Entrez Gene and returns them as a list.
        '''
        # Limited to IDs which are current and not obsolete (i.e. "alive")
        term = '("alive"[Properties]) AND %s[Taxonomy ID]' % host

        E = Entrez.esearch(db="gene", term=term, retmode="text",
                           retmax=1000000)
        res = Entrez.read(E)
        E.close()

        return res['IdList']


class EntrezTaxonomyAnnotation(EntrezAnnotation):
    def __init__(self, outdb, email):
        EntrezAnnotation.__init__(self, "entreztaxonomy", outdb, email)

    def download(self, ids):
        E = Entrez.efetch(db="taxonomy", id=ids, retmode="xml", retmax=1000000)
        res = Entrez.read(E)
        self.dataset = res


class MyGeneInfoAnnotation(APIAnnotation):
    def __init__(self, prefix, source, outdb, options=None, ohost=None,
                 email=None):
        APIAnnotation.__init__(self, prefix, source, outdb, options, ohost,
                               email)

    def test(self):
        '''
        Test MyGeneInfo API - Look up symbol for Entrez ID 60489, result
        should be APOBEC3G.
        '''
        h = httplib2.Http()
        res, con = h.request(self.datasource, 'POST',
                             'ids=60489&fields=symbol',
                             headers={'content-type':
                                      'application/x-www-form-urlencoded'})
        # convert json format into Python
        jcon = json.loads(con)
        if jcon[0]['symbol'] == "APOBEC3G":
            return 1
        else:
            return 0

    def download(self, idlist, fields):
        '''
        Gets information from MyGeneInfo.org for the list of genes pulled from
        the Entrez Gene database.  This information is combined into a big
        python dictionary.
        The fields to get from MyGeneInfo are specified in pipeline.ini, the
        list of all possible fields is at
        http://docs.mygene.info/en/latest/doc/data.html#available-fields
        Each MyGeneInfo 'object' is formatted slightly differently, depending
        on the field and host, so adding another field means writing a new
        function to parse the information retrieved (similar to the
        ParseAndLoad functions below).
        '''

        DB = dict()
        # If the list of entrez ids is long, cut it into chunks of 500 ids to
        # send to the API in each query
        # (this size seems to run the fastest overall).
        if len(idlist) >= 500:
            ids = np.array(idlist)
            a = len(ids) / 500
            idsets = np.array_split(ids, a)
        else:
            idsets = [idlist]
        for idset in idsets:
            h = httplib2.Http()
            headers = {'content-type': 'application/x-www-form-urlencoded'}

            # get all the ids and fields specified
            params = 'ids=%s&fields=%s' % (",".join(idset), ",".join(fields))
            res, con = h.request(self.datasource, 'POST',
                                 params, headers=headers)

            # convert JSON format to nested Python dictionaries.
            jcon = json.loads(con)
            for line in jcon:
                for f in line.keys():
                    if f not in DB:
                        DB[f] = dict()
                    DB[f][line['query']] = line[f]
        self.dataset = DB


class SymbolAnnotation(MyGeneInfoAnnotation):
    def __init__(self, source, outdb, host, shost):
        MyGeneInfoAnnotation.__init__(self, "symbol", source, outdb,
                                      ohost=host)
        self.shost = shost

    def parse(self):
        entrez, symbol = [], []
        res = self.dataset['symbol']
        for item in res.items():
            entrez.append(item[0])
            symbol.append(item[1])
        df = pd.DataFrame(zip(entrez, symbol), columns=['entrez', 'symbol'])
        df = self.translate(df, 'entrez', 'ensemblg')
        self.resultspd['ensemblg2symbol_%s$geneid' % self.shost] = df
        self.storeTranslation(df['ensemblg'], df['symbol'],
                              ['ensemblg', 'symbol'],
                              "ensemblg2symbol_%s.tsv" % self.shost)


class EnsemblAnnotation(MyGeneInfoAnnotation):
    def __init__(self, source, outdb):
        MyGeneInfoAnnotation.__init__(self, "ensembl", source, outdb)

    def parse(self):
        ensembldict = self.dataset['ensembl']

        # store the relationships between various gene ids as lists of tuples
        entrez2ensembl = []
        ensembl2entrez = []
        g2t = []
        t2g = []
        g2p = []
        p2g = []

        # for each entrez ID in the dictionary of ensembl results from
        # mygeneinfo
        for q in ensembldict:
            # each ID leads to another dictionary with keys 'gene',
            # 'transcript'
            # and 'protein' containing ensembl ids of these types.
            res = ensembldict[q]
            # in the "ensembldict" created from the mygeneinfo output
            # 1:1 mappings have dictionary values stored as strings
            # and 1 to many mappings have dictionary values stored as
            # lists.  This makes a list with one element if the value is a
            # string, so that subsequent steps can be applied in both cases.
            if type(res) is not list:
                res = [res]
            for L in res:
                # L['gene'], L['transcript'] and L['protein'] are also a
                # mixture of lists and strings depending if there is one or
                # more than one value, most of the following is to deal
                # with this
                if type(L['gene']) is list:
                    entrez2ensembl += [q] * len(L['gene'])
                    ensembl2entrez += L['gene']
                    glist = L['gene']
                else:
                    entrez2ensembl.append(q)
                    ensembl2entrez.append(L['gene'])
                    glist = [L['gene']]
                for gene in glist:
                    if 'transcript' in L:
                        if type(L['transcript']) is list:
                            g2t += [gene] * len(L['transcript'])
                            t2g += L['transcript']
                        else:
                            g2t.append(gene)
                            t2g.append(L['transcript'])
                    if 'protein' in L:
                        if type(L['protein']) is list:
                            g2p += [gene] * len(L['protein'])
                            p2g += L['protein']
                        else:
                            g2p.append(gene)
                            p2g.append(L['protein'])
        self.resultsz['ensemblg2entrez$geneid'] = [zip(entrez2ensembl,
                                                       ensembl2entrez),
                                                   ['ensemblg', 'entrez'],
                                                   ['int', 'varchar(25)']]

        self.resultsz['ensemblg2ensemblt$other'] = [zip(g2t, t2g),
                                                    ['ensemblg', 'ensemblt'],
                                                    ['varchar(25)',
                                                     'varchar(25)']]
        self.resultsz['ensemblg2ensemblp$other'] = [zip(g2p, p2g),
                                                    ['ensemblg', 'ensemblp'],
                                                    ['varchar(25)',
                                                     'varchar(25)']]
        self.storeTranslation(ensembl2entrez, entrez2ensembl, ['ensemblg',
                                                               'entrez'],
                              'ensemblg2entrez.tsv')


class GoAnnotation(MyGeneInfoAnnotation):
    def __init__(self, source, outdb, options):
        MyGeneInfoAnnotation.__init__(self, "go", source, outdb)
        self.options = options

    def parse(self):
        gores = self.dataset['go']
        godict = dict()
        # options from pipeline.ini for which go namespaces are of interest
        options = self.options.split(",")
        goids = []
        goset = set()
        allgoids = []
        entrez2go = []
        gotypes = []
        # for each entrez ID in the go dictionary
        for q in gores:
            # res is a dictionary with keys BP, CC and MF
            res = gores[q]
            for cat in ['BP', 'MF', 'CC']:
                if cat in res and (cat in options or 'all' in options):
                    # As for ensembl above, this deals with the mixture of
                    # lists and strings as dictionary values returned by
                    # mygeneinfo
                    if type(res[cat]) is not list:
                        res[cat] = [res[cat]]
                    for L in res[cat]:
                        allgoids.append(L['id'])
                        entrez2go.append(q)
                        if L['id'] not in goset:
                            goset = goset | set([L['id']])
                            goids.append(L['id'])
                            gotypes.append(cat)
                            for subterm in ['term', 'evidence']:
                                if subterm not in godict:
                                    godict[subterm] = []
                                if subterm in L:
                                    godict[subterm].append(L[subterm])
                                else:
                                    godict[subterm].append("")
        dfs = pd.DataFrame(zip(entrez2go, allgoids), columns=['entrez', 'go'])
        results = self.translate(dfs, 'entrez', 'ensemblg')
        self.resultspd['ensemblg2go$annot'] = results
        self.resultsz['go$details'] = [zip(goids, gotypes, godict['term'],
                                           godict['evidence']),
                                       ['go',
                                        'gotype',
                                        'goterm',
                                        'goevidence'],
                                       ['varchar(25)',
                                        'varchar(25)',
                                        'varchar(255)',
                                        'varchar(25)']]


class PathwayAnnotation(MyGeneInfoAnnotation):
    def __init__(self, source, outdb, options):
        MyGeneInfoAnnotation.__init__(self, "pathway", source, outdb)
        self.options = options

    def parse(self):
        '''
        Parses the "pathway" dictionary extracted from MyGeneInfo above.
        The pathway annotation sets of interest are specified in pipeline.ini,
        or if 'all' is specified all pathway annotations available in
        the database are used.
        Each annotation set is used to generate 2 database tables:
        entrez2xxx - translates entrez ids to annotation ids
        xxx - annotation ids and the pathways they correspond to
        '''
        DB = self.dataset
        if 'pathway' not in DB:
            print "pathway annotations not available for this species"
            return 0
        pathwaydict = DB['pathway']
        D = dict()

        # If specific annotation sets to use are provided in pipeline.ini
        # use those
        # otherwise get all the keys from the dictionaries in pathwaydict.
        typelist = self.options.split(",")

        if 'all' in typelist:
            allkeys = set()
            for k in DB['pathway'].keys():
                allkeys = allkeys | set(DB['pathway'][k].keys())
            typelist = list(allkeys)

        # Make empty lists to store the lists of tuples to load into the
        # database
        # The D dictionary contains four lists per host, named host_ids2entrez,
        # host_entrez2ids, host_id and host_name.
        for db in typelist:
            D["%s_ids2entrez" % db] = []
            D["%s_entrez2ids" % db] = []
            D["%s_id" % db] = []
            D["%s_name" % db] = []

        # for each entrez id in the pathway dictionary
        for q in pathwaydict.keys():
            # for each of the requested annotation types
            for db in typelist:
                idset = set()
                if db in pathwaydict[q]:
                    results = pathwaydict[q][db]
                    if type(results) is not list:
                        results = [results]
                    for b in results:
                        D['%s_ids2entrez' % db].append(b['id'])
                        D["%s_entrez2ids" % db] .append(q)
                        if b['id'] not in idset:
                            D["%s_id" % db].append(b['id'])
                            D["%s_name" % db].append(b['name'])
                            idset.add(b['id'])

        # put these into the database
        for db in typelist:
            ids2entrez = D["%s_ids2entrez" % db]
            entrez2ids = D["%s_entrez2ids" % db]
            id = D["%s_id" % db]
            name = D["%s_name" % db]
            identrez = zip(entrez2ids, ids2entrez)
            idname = zip(id, name)
            dfs = pd.DataFrame(identrez, columns=['entrez', db])
            self.translate(dfs, 'entrez', 'ensemblg')
            self.resultspd['ensemblg2%s$annot' % db] = dfs
            self.resultsz['%s$details' % db] = [idname, [db, 'term'],
                                                ['varchar(25)',
                                                 'varchar(25)']]


class HomologeneAnnotation(MyGeneInfoAnnotation):
    def __init__(self, source, outdb, options, ohost, email):
        MyGeneInfoAnnotation.__init__(self, "homologene", source, outdb,
                                      options, ohost, email)

    def parse(self):
        '''
        Uses information from the homologene dataset retrieved from MyGeneInfo
        to link entrez IDs from the host organism to those of homologous genes
        in other organisms of interest.
        The organisms of interest can be specified in pipeline.ini or, if 'all'
        is specified, every available organism is used (~25 organisms).
        A summary table "species_ids" is generated showing which
        host taxonomy id corresponds to which host scientific name.
        For each taxon a table is generated, taxonname_entrez, linking host
        entrez ids to those of the taxon.
        Entrez ids for each taxon are translated into gene symbols, these are
        stored in tables as entrez_taxonname_symbol.
        '''
        homodb = self.dataset['homologene']
        taxa = str(self.options).split(",")
        if '9606' not in taxa:
            taxa.append('9606')

        # retrieve IDs for all the taxa listed in
        # PARAMS["my_gene_info_homologene"]
        # or if this is "all" check which taxa are available

        if 'all' in taxa:
            allspp = set()
            for q in homodb.keys():
                geneset = homodb[q]
                if type(geneset) is not list:
                    geneset = [geneset]
                for L in geneset:
                    spp = [str(line[0]) for line in L['genes']]
                    allspp = allspp | set(spp)
            allspp = list(allspp)
        else:
            taxa = set(taxa)
            taxa.add(str(self.ohost))
            allspp = list(taxa)
        # retrieve species names of the taxa from entrez taxonomy

        asp = ",".join(allspp)
        E = EntrezTaxonomyAnnotation(self.outdb, self.email)
        E.download(asp)
        res = E.dataset

        names = dict()
        for r in res:
            names[r['TaxId']] = "_".join(r['ScientificName'].split(" "))

        entrezdict = dict()
        for spp in allspp:
            entrezdict[spp] = (([], []))

        i = 0
        for q in homodb.keys():
            for item in homodb[q]['genes']:
                hostid = str(item[0])
                if hostid in allspp:
                    hostgeneid = str(item[1])
                    entrezdict[hostid][0].append(q)
                    entrezdict[hostid][1].append(hostgeneid)
            i += 1
        for host in entrezdict:
            R = zip(entrezdict[host][0], entrezdict[host][1])
            R = pd.DataFrame(R, columns=['entrez', 'entrez_%s' % host])
            R = self.translate(R, 'entrez', 'ensemblg')
            Sym = SymbolAnnotation('http://mygene.info/v2/gene', self.outdb,
                                   host=host, shost=names[host])
            Sym.download(entrezdict[host][1], ['symbol'])
            Sym.parse()
            df = Sym.dataset['symbol']
            df = [((key, df[key])) for key in df]
            df = pd.DataFrame(df, columns=['entrez_%s' % host,
                                           'symbol_%s' % host])
            R['entrez_%s' % host] = R['entrez_%s' % host].astype(int)
            df['entrez_%s' % host] = df['entrez_%s' % host].astype(int)
            df = df.merge(R,
                          left_on='entrez_%s' % host,
                          right_on='entrez_%s' % host)
            df = df.drop('entrez_%s' % host, 1)
            self.resultspd['ensemblg2symbol_%s$geneid' % names[host]] = df
            self.storeTranslation(df['ensemblg'],
                                  df['symbol_%s' % host],
                                  ['ensemblg', 'symbol_%s' % names[host]],
                                  'ensemblg2symbol_%s.tsv' % names[host])


class DataMineAnnotation(APIAnnotation):
    def __init__(self, prefix, source, outdb, chost, ind, views, constraints,
                 hostid):
        APIAnnotation.__init__(self, prefix, source, outdb)
        self.chost = chost
        self.ind = ind
        self.views = views
        self.hostid = hostid
        self.constraints = constraints

    def test(self):
        service = Service('http://www.humanmine.org/humanmine/service')
        query = service.new_query("Gene")
        query.add_view("symbol")
        query.add_constraint("Gene", "LOOKUP", "APOBEC3G", code="A")
        for row in query.rows():
            symbol = row['symbol']
        if symbol == "APOBEC3G":
            return 1
        else:
            return 0

    def download(self, genes, fields):
        constraints = self.constraints
        views = self.views
        glist = np.array(genes)
        if len(glist) > 1000:
            a = len(glist) / 1000
            segs = np.array_split(glist, a)
        else:
            segs = [glist]

        # store the data in here
        z = []

        # API uses letters to distinguish between constraints
        alpha = list(string.ascii_uppercase)

        for seg in segs:
            # Connect to the API
            service = Service(self.datasource)
            query = service.new_query("Gene")
            query.add_view(",".join(views))
            # Some databases require a host name
            if self.hostid != "":
                query.add_constraint("Gene", "LOOKUP", ",".join(seg),
                                     self.hostid, code="A")
            else:
                query.add_constraint("Gene", "LOOKUP", ",".join(seg), code="A")

            # Apply the constraints
            if len(constraints) != 0:
                i = 1
                for constraint in constraints:
                    letter = alpha[i]
                    if len(constraint.split("=")) == 2:
                        L = constraint.split("=")
                        query.add_constraint(L[0], "=", L[1], code=letter)
                    elif re.search("IS NOT NULL", constraint):
                        p1 = constraint.replace(" IS NOT NULL", "")
                        query.add_constraint(p1, "IS NOT NULL", code=letter)
                    i = i + 1

            # Parse the output into a list of tuples
            j = 0
            for row in query.rows():
                t = [row['symbol']]
                for v in views:
                    t.append(row[v])
                z.append(tuple(t))
                j += 1
        self.dataset = z

    def parse(self):
        z = self.dataset
        symbols = []
        terms = []
        other = []
        ind = self.ind
        for k in z:
            symbols.append(k[0])
            terms.append(k[ind+1])
            L = list(k)
            L.remove(k[0])
            L.remove(k[ind+1])
            k2 = [k[ind+1]]
            k2 += L
            other.append(tuple(k2))

        dfs = pd.DataFrame(zip(symbols, terms),
                           columns=['symbol_%s' % self.chost, 'term'])
        dfs = self.translate(dfs, 'symbol_%s' % self.chost, 'ensemblg')
        views = self.views
        cols = [v.replace(".", "_") for v in views[0:ind] + views[(ind+1):]]
        cols = ['term'] + cols
        other = pd.DataFrame(other, columns=cols)
        self.resultspd["%s$details" % self.prefix] = other
        self.resultspd["ensemblg2%s$annot" % self.prefix] = dfs


class MGIAnnotation(DataMineAnnotation):
    def __init__(self, source, outdb):
        prefix = "mgi"
        chost = "Mus_musculus"
        ind = 3
        views = ["ontologyAnnotations.ontologyTerm.description",
                 "ontologyAnnotations.ontologyTerm.name",
                 "ontologyAnnotations.ontologyTerm.namespace",
                 "ontologyAnnotations.ontologyTerm.identifier"]
        constraints = [
            ("Gene.ontologyAnnotations.ontologyTerm.namespace=MPheno.ontology")]
        hostid = 'M. musculus'
        DataMineAnnotation.__init__(self, prefix, source, outdb,
                                    chost, ind, views, constraints, hostid)


class MousePathwayAnnotation(DataMineAnnotation):
    def __init__(self, source, outdb):
        prefix = "mousepathway"
        chost = "Mus_musculus"
        ind = 0
        views = ["pathways.identifier", "pathways.name"]
        constraints = []
        hostid = 'M. musculus'
        DataMineAnnotation.__init__(self, prefix, source, outdb,
                                    chost, ind, views, constraints, hostid)


class HPOAnnotation(DataMineAnnotation):
    def __init__(self, source, outdb):
        prefix = "hpo"
        chost = "Homo_sapiens"
        ind = 3
        views = ["diseases.hpoAnnotations.hpoTerm.description",
                 "diseases.hpoAnnotations.hpoTerm.name",
                 "diseases.hpoAnnotations.hpoTerm.namespace",
                 "diseases.hpoAnnotations.hpoTerm.identifier"]
        constraints = []
        hostid = ''
        DataMineAnnotation.__init__(self, prefix, source, outdb,
                                    chost, ind, views, constraints, hostid)


class OntologyAnnotation(APIAnnotation):
    '''
    Parses an OWL format ontology to find the hierarchy of terms.
    Generates a database table, tablename with a list of terms in column
    1 and all the terms of which they are a direct subclass in column 2, with
    one term per row.
    For example, from the human phenotype ontology:
    Glomerulonephritis (HP:0000099) is a subclass of nephritis (HP:0000123)
    and of "Abnormality of the Glomerulus" (HP:0000095).
    This would be recorded as:
    id    subclass_of
    HP:0000099    HP:0000095
    HP:0000099    HP:0000123
    '''
    def __init__(self, prefix, datasource, outdb):
        APIAnnotation.__init__(self, prefix, datasource, outdb)

    def test(self):
        try:
            urllib2.urlopen('http://purl.obolibrary.org/obo/go.owl',
                            timeout=1)
            return True
        except urllib2.URLError:
            pass

    def download(self, genes=None, fields=None):
        # download an up to date ontology file, parse the xml data into a
        # Python "ElementTree" and delete the ontology file.
        ontologyfile = P.getTempFilename(".")
        os.system("wget -O %s %s" % (ontologyfile, self.datasource))
        tree = ET.parse(ontologyfile)
        os.remove(ontologyfile)
        self.dataset = tree

    def parse(self):
        # Traverse the tree and identify the classes and subclasses
        # The ns dictionary links "prefixes" of xml tags to the full tag
        # names, see here for more info -
        # https://docs.python.org/2/library/xml.etree.elementtree.html
        ns = {'owl': 'http://www.w3.org/2002/07/owl#',
              'rdfs': 'http://www.w3.org/2000/01/rdf-schema#'}
        tree = self.dataset
        root = tree.getroot()
        r = dict()
        for aclass in root.findall('owl:Class', ns):
            nam = aclass.attrib.values()[0].split("/")[-1]
            for subclassof in aclass.findall('rdfs:subClassOf', ns):
                if nam not in r:
                    r[nam] = set()
                if len(subclassof.attrib.values()) == 1:
                    r[nam].add(subclassof.attrib.values()[0].split("/")[-1])

        klist = []
        vlist = []
        for item in r.items():
            key = item[0]
            vals = item[1]
            for v in vals:
                klist.append(key.replace("_", ":"))
                vlist.append(v.replace("_", ":"))
        self.resultsz['%s$ont' % self.prefix] = [zip(klist, vlist),
                                                 ['id', 'subclass_of'],
                                                 ['varchar(250)',
                                                  'varchar(250)']]


@cluster_runnable
def runall(obj, genelist=None, fields=None):
    obj.runall(genelist, fields)










