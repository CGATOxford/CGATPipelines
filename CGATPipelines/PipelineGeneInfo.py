import CGATPipelines.Pipeline as P
from Bio import Entrez
import numpy as np
import httplib2
import json as json
import sqlite3
from intermine.webservice import Service as SS
import string
import re
import os
import xml.etree.ElementTree as ET
import pandas as pd
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from future.moves.urllib.request import urlopen
from CGATPipelines.Pipeline import cluster_runnable
import mygene


def readGeneList(filename):
    '''
    Reads a file containing one gene ID per line into a list
    '''
    return [line.strip() for line in IOTools.openFile(filename).readlines()]


def getSymbols(hostfile):
    '''
    Reads a file containing pairs of gene symbols and another ID type and
    returns a list of gene symbols
    '''
    return [line.strip().split("\t")[0]
            for line in IOTools.openFile(hostfile).readlines()]


class APIAnnotation(object):
    '''
    Parent class for any annotation type to be downloaded via an API,
    parsed and loaded into a database.
    '''

    def __init__(self, prefix, datasource, outdb, options=None, ohost=None,
                 email=None):
        '''
        prefix - prefix used to recognise which annotation is being performed
        datasource - link used for API
        outdb - output database
        options - optional string to be parsed by the APIAnnotation subclass
        ohost - original species of interest used by some subclasses
        email - entrez requires an email address to use its API
        '''
        self.prefix = prefix
        self.datasource = datasource
        self.outdb = outdb
        self.options = options
        self.ohost = ohost
        self.email = email
        self.resultsz = dict()
        self.resultspd = dict()

    def runall(self, genelist, fields, scope, species):
        '''
        Basic method to take the details of an API annotation and:
        - test connectivity to the API
        - download the data
        - parse the data into a useable format
        - load the data into a database
        '''
        T = self.test()
        if T != 1:
            E.warn("API test failed for %s" % self.prefix)
            return None

        self.download(genelist, fields=fields, scope=scope, species=species)
        self.parse()
        self.loadPDTables()
        self.loadZippedTables()

    def translate(self, dataframe, tfrom, tto):
        '''
        Translates a column of geneids in a dataframe from type
        tfrom to type tto
        Requires both ID types to have been downloaded in a previous
        step - uses an existing stored dataframe with tfrom and tto columns
        '''
        assert ((os.path.exists('%s2%s.tsv' % (tfrom, tto))) |
                (os.path.exists('%s2%s.tsv' % (tto, tfrom)))), """
                translation from %s to %s needs to be available prior
                to %s parsing""" % (tfrom, tto, self.prefix)
        try:
            ensdf = self.getTranslation('%s2%s.tsv' % (tfrom, tto))
        except:
            ensdf = self.getTranslation('%s2%s.tsv' % (tto, tfrom))
        # dataframe is the pandas dataframe to be translated
        dataframe[tfrom] = dataframe[tfrom].astype(str)
        # ensdf is the pandas dataframe used to translate
        ensdf[tfrom] = ensdf[tfrom].astype(str)
        # use pandas merge to append the columns from ensdf to dataframe
        # then remove unneeded columns
        return dataframe.merge(
            ensdf, left_on=tfrom, right_on=tfrom).drop(tfrom, 1)

    def loadZippedTables(self):
        '''
        Adds a table to the database for each key in the self.resultsz
        dictionary.
        Each key is the name of the database table and consists of a list of
        three items:
        self.resultsz[key][0] - a list of tuples containing the data to load,
        where each tuple corresponds to one row and each item within the tuple
        one column
        self.resultsz[key][1] - table column names as a Python list
        self.resultsz[key][2] -  column types as strings in a Python list

        E.g.
        self.resultsz[key] = [[(1111, 'aaa'), (2222, 'bbb'), (3333, 'ccc)],
                              ["entrez_id", "gene_symbol"],
                              ['int', 'varchar (25)']]

        would create the table "entrez2symbol" in the database with the
        layout
        entrez_id    gene_symbol
        1111         aaa
        2222         bbb
        3333         ccc


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
            statement = 'CREATE TABLE IF NOT EXISTS %s (%s);' % (
                tablename, cnamestypes)
            c.execute(statement)

            # populate the table with the data
            statement = '''INSERT INTO %s (%s) VALUES (%s);''' % (tablename,
                                                                  cnamestring,
                                                                  qs)
            outf.write("%s\n" % statement)

            c.executemany(statement, zipped)
            outf.write("%s\n" % statement)
            dbh.commit()
            outf.close()
        # disconnect
        c.close()
        dbh.close()

    def loadPDTables(self):
        '''
        Loads a table to the database corresponding to every key in the
        self.resultspd dictionary.
        Keys are the names of the tables and values are the data to upload.
        '''
        dbh = sqlite3.connect(self.outdb)
        for tabkey in self.resultspd:
            tab = self.resultspd[tabkey]
            tab.to_sql(tabkey, dbh, index=False, if_exists='replace')
            tab.to_csv("%s.load" % tabkey, sep="\t", index=False,
                       encoding='utf-8')
        dbh.close()

    def storeTranslation(self, genelist1, genelist2, cnames, out):
        '''
        When an object is parsed containing a translation from one type
        of gene id to another (e.g. entrez ID to ensemblg ID) it is sometimes
        useful to have it available easily later - this stores a basic
        tab delimited file containing this information.
        '''
        df = pd.DataFrame(list(zip(genelist1, genelist2)), columns=cnames)
        df.to_csv(out, sep="\t", index=False)

    def getTranslation(self, translation):
        '''
        Reads the tab delimited files generated by storeTranslation and
        regenerates the pandas dataframe.
        '''
        df = pd.read_csv(translation, sep="\t")
        return df


class EntrezAnnotation(APIAnnotation):
    '''
    Used to read, parse and load data from the Entrez API, which is accessed
    via the Entrez BioPython package.
    '''

    def __init__(self, prefix, outdb, email):
        '''
        Prefix - used to reconise the annotation currently running
        outdb - output database
        email - entrez requires an email address to use its API, this
        should be provided in the pipeline.ini file.
        '''
        APIAnnotation.__init__(self, prefix, "Entrez", outdb)
        Entrez.email = email

    def test(self):
        '''
        Test Entrez API is connecting and working.
        Looks up symbol APOBEC3G, Entrez ID 60489 should be
        amongst the results
            '''
        Ent = Entrez.esearch(db="gene", term="(APOBEC3G[Preferred+Symbol])",
                             retmode="text", retmax=1000000)
        res = Entrez.read(Ent)
        Ent.close()
        if "60489" in res['IdList']:
            return 1
        else:
            return 0


class EntrezGeneAnnotation(EntrezAnnotation):
    '''
    Entrez Gene annotation object
    '''

    def __init__(self, outdb, email):
        EntrezAnnotation.__init__(self, "entrezgene", outdb, email)

    def download_all(self, host, count=1000000000):
        '''
        Gets all the Gene IDs for a particular host, specified in PARAMS, from
        Entrez Gene and returns them as a list.
        '''
        # Limited to IDs which are current and not obsolete (i.e. "alive")
        term = '("alive"[Properties]) AND %s[Taxonomy ID]' % host

        Ent = Entrez.esearch(db="gene", term=term, retmode="text",
                             retmax=count)
        res = Entrez.read(Ent)
        Ent.close()

        return res['IdList']


class EntrezTaxonomyAnnotation(EntrezAnnotation):
    '''
    Entrez Taxonomy Annotation object
    '''

    def __init__(self, outdb, email):
        EntrezAnnotation.__init__(self, "entreztaxonomy", outdb, email)

    def download(self, ids):
        '''
        Fetch data from the Entrez Taxonomy database for the Taxonomy IDs in
        ids
        '''
        Ent = Entrez.efetch(db="taxonomy",
                            id=ids, retmode="xml", retmax=1000000)
        res = Entrez.read(Ent)
        self.dataset = res


class MyGeneInfoAnnotation(APIAnnotation):
    '''
    Used to download, parse and load data from the mygene.info API.
    Data is provided in nested dictionaries, the number of levels
    and data types are not consistent for different annotations so
    a different parser is currently needed for each annotation type.
    '''

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
        mg = mygene.MyGeneInfo()

        con = mg.querymany(['60489'], scope='entrezgene', fields='symbol')

        if con[0]['symbol'] == "APOBEC3G":
            return 1
        else:
            return 0

    def download(self, idlist, fields, scope, species):
        '''
        Gets information from MyGeneInfo.org for the list of genes pulled from
        the Entrez Gene database.  This information is combined into a big
        python dictionary.
        The fields to get from MyGeneInfo are specified in pipeline.ini, the
        list of all possible fields is at
        http://docs.mygene.info/en/latest/doc/data.html#available-fields
        Each MyGeneInfo 'object' is formatted slightly differently, depending
        on the field and host, so adding another field means writing a new
        function to parse the information retrieved.
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

        mg = mygene.MyGeneInfo()
        for idset in idsets:
            con = mg.querymany(idset, scope=scope, fields=fields,
                               species=species)
            for line in con:
                for f in list(line.keys()):
                    if f not in DB:
                        DB[f] = dict()
                    if self.prefix == 'ensembl':
                        DB[f].setdefault(line['query'], [])
                        DB[f][line['query']].append(line[f])
                    else:
                        DB[f][line['query']] = line[f]
        self.dataset = DB


class SymbolAnnotation(MyGeneInfoAnnotation):
    '''
    Gene symbol annotation from mygene.info
    '''

    def __init__(self, source, outdb, host, shost):
        '''
        host is the entrez taxonomy ID of the species to download symbols for
        shost is the scientific name of this species, seperated by "_" e.g.
        Homo_sapiens
        '''
        MyGeneInfoAnnotation.__init__(self, "symbol", source, outdb,
                                      ohost=host)
        self.shost = shost

    def parse(self):
        '''
        Parse the symbol annotation stored by self.download as
        self.dataset['symbol'] into a pandas dataframe with one line
        per unique gene ID (either ensemblg or entrez).
        Store the translation as a tsv formatted file to use later.
        '''
        entrez, symbol = [], []
        res = self.dataset['symbol']
        for item in list(res.items()):
            entrez.append(item[0])
            symbol.append(item[1])
        df = pd.DataFrame(list(zip(entrez, symbol)),
                          columns=['entrez', 'symbol_%s' % self.ohost])
        # df = self.translate(df, 'entrez', 'ensemblg')
        self.resultspd['entrez2symbol_%s$geneid' % self.ohost] = df
        self.storeTranslation(df['entrez'], df['symbol_%s' % self.ohost],
                              ['entrez', 'symbol_%s' % self.ohost],
                              "entrez2symbol_%s.tsv" % self.ohost)


class EnsemblAnnotation(MyGeneInfoAnnotation):
    '''
    Used to parse Ensembl data from mygene.info
    '''

    def __init__(self, source, outdb, ohost):
        MyGeneInfoAnnotation.__init__(self, "ensembl", source, outdb,
                                      ohost=ohost)

    def parse(self):
        '''
        Parse the dictionary of dictionaries stored by self.download() in
        self.dataset['ensembl'] into lists of tuples to upload into the
        database
        '''
        ensembldict = self.dataset['ensembl']

        # store the relationships between various gene ids as lists of tuples
        sym2ensembl = []
        ensembl2sym = []
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
            resL = ensembldict[q]
            # in the "ensembldict" created from the mygeneinfo output
            # 1:1 mappings have dictionary values stored as strings
            # and 1 to many mappings have dictionary values stored as
            # lists.  This makes a list with one element if the value is a
            # string, so that subsequent steps can be applied in both cases.
            for res in resL:
                if type(res) is not list:
                    res = [res]
                for L in res:
                    # L['gene'], L['transcript'] and L['protein'] are also a
                    # mixture of lists and strings depending if there is one or
                    # more than one value, most of the following is to deal
                    # with this
                    if type(L['gene']) is list:
                        sym2ensembl += [q] * len(L['gene'])
                        ensembl2sym += L['gene']
                        glist = L['gene']
                    else:
                        sym2ensembl.append(q)
                        ensembl2sym.append(L['gene'])
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
        # load everything into the self.resultsz dictionary so that
        # self.loadZippedTables will put it into the database
        self.resultsz['ensemblg2symbol_%s$geneid' % self.ohost] = [list(
            zip(ensembl2sym,
                sym2ensembl)), ['ensemblg',
                                'symbol_%s' % self.ohost],
                                ['int', 'varchar(25)']]

        self.resultsz['ensemblg2ensemblt$other'] = [list(zip(g2t, t2g)),
                                                    ['ensemblg', 'ensemblt'],
                                                    ['varchar(25)',
                                                     'varchar(25)']]
        self.resultsz['ensemblg2ensemblp$other'] = [list(zip(g2p, p2g)),
                                                    ['ensemblg', 'ensemblp'],
                                                    ['varchar(25)',
                                                     'varchar(25)']]
        # store the entrez to ensemblg translation for later
        self.storeTranslation(ensembl2sym,
                              sym2ensembl, ['ensemblg',
                                            'symbol_%s' % self.ohost],
                              'ensemblg2symbol_%s.tsv' % self.ohost)
        ens2symdf = pd.DataFrame(zip(ensembl2sym, sym2ensembl),
                                 columns=['ensemblg',
                                          'symbol_%s' % self.ohost])
        res = self.translate(ens2symdf, 'symbol_%s' % self.ohost,
                             "entrez")
        self.storeTranslation(res['ensemblg'], res['entrez'],
                              ['ensemblg', 'entrez'], "ensemblg2entrez.tsv")

        self.resultsz['ensemblg2entrez$geneid'] = [list(zip(res['ensemblg'],
                                                            res['entrez'])),
                                                   ['ensemblg', 'entrez'],
                                                   ['varchar(25)',
                                                    'varchar(25)']]


class GoAnnotation(MyGeneInfoAnnotation):
    '''
    GO gene ontology annotation downloaded from mygene.info
    '''

    def __init__(self, source, outdb, options, ohost):
        '''
        options can be BP, CC, MF or any combination of these, stored as a
        string and comma delimited (from pipeline.ini)
        '''
        MyGeneInfoAnnotation.__init__(self, "go", source, outdb, ohost=ohost)
        self.options = options

    def parse(self):
        '''
        Parses the dictionary of dictionaries stored by self.download() into
        lists of tuples to load into the database.
        '''
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
        dfs = pd.DataFrame(list(zip(entrez2go, allgoids)),
                           columns=['symbol_%s' % self.ohost, 'go'])
        results = self.translate(dfs, 'symbol_%s' % self.ohost, 'ensemblg')
        # add the ensemblg to go term dataframe to resultspd for
        # self.loadPDTables to load into the database
        self.resultspd['ensemblg2go$annot'] = results
        # add the metadata to the resultsz dictionary for
        # self.loadZippedTables() to load into the database.
        self.resultsz['go$details'] = [list(zip(goids, gotypes, godict['term'],
                                                godict['evidence'])),
                                       ['go',
                                        'gotype',
                                        'goterm',
                                        'goevidence'],
                                       ['varchar(25)',
                                        'varchar(25)',
                                        'varchar(255)',
                                        'varchar(25)']]


class PathwayAnnotation(MyGeneInfoAnnotation):
    '''
    Pathway Annotations from mygene.info.  The 'pathway' dictionary returned
    by self.parse() contains annotations from a number of pathway
    database, those of interest (or
    all) can be specified in the pipeline.ini.
    '''

    def __init__(self, source, outdb, options, ohost):
        '''
        options - which pathway databases to keep (as a comma delimited
        string from pipeline.ini)
        '''
        MyGeneInfoAnnotation.__init__(self, "pathway", source, outdb,
                                      ohost=ohost)
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
            print("pathway annotations not available for this species")
            return 0
        pathwaydict = DB['pathway']
        D = dict()

        # If specific annotation sets to use are provided in pipeline.ini
        # use those
        # otherwise get all the keys from the dictionaries in pathwaydict.
        typelist = self.options.split(",")

        if 'all' in typelist:
            allkeys = set()
            for k in list(DB['pathway'].keys()):
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
        for q in list(pathwaydict.keys()):
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

        # put these into self.resultsz and self.resultspd for
        # self.loadxxx to put into the database
        for db in typelist:
            ids2entrez = D["%s_ids2entrez" % db]
            entrez2ids = D["%s_entrez2ids" % db]
            id = D["%s_id" % db]
            name = D["%s_name" % db]
            identrez = list(zip(entrez2ids, ids2entrez))
            idname = list(zip(id, name))
            dfs = pd.DataFrame(identrez, columns=['symbol_%s' % self.ohost,
                                                  db])
            dfs = self.translate(dfs, 'symbol_%s' % self.ohost, 'ensemblg')
            self.resultspd['ensemblg2%s$annot' % db] = dfs
            self.resultsz['%s$details' % db] = [idname, [db, 'term'],
                                                ['varchar(25)',
                                                 'varchar(25)']]


class HomologeneAnnotation(MyGeneInfoAnnotation):
    '''
    Uses information from the homologene dataset retrieved from MyGeneInfo
    to link entrez IDs from the host organism to those of homologous genes
    in other organisms of interest.
    The organisms of interest can be specified in pipeline.ini or, if 'all'
    is specified, every available organism is used ($25 organisms).
    '''

    def __init__(self, source, outdb, options, ohost, email):
        MyGeneInfoAnnotation.__init__(self, "homologene", source, outdb,
                                      options, ohost, email)

    def parse(self):
        '''
        Entrez ids for each taxon are translated into gene symbols, these are
        stored in tables as ensemblg2taxonname_symbol.
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
            for q in list(homodb.keys()):
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
        Ent = EntrezTaxonomyAnnotation(self.outdb, self.email)
        Ent.download(asp)
        res = Ent.dataset

        # translate taxonomy IDs into scientific names
        names = dict()
        for r in res:
            names[r['TaxId']] = "_".join(r['ScientificName'].split(" "))

        # entrezdict will contain two lists of each species
        # list 0 is the entrez ID in the original species
        # list 1 is the entrez ID in the species with the homologue
        entrezdict = dict()
        for spp in allspp:
            entrezdict[spp] = (([], []))

        i = 0
        # for each entrez id in the original species
        for q in list(homodb.keys()):
            # retrieve the list of genes (stored as a list of lists as
            # [[host, id], [host, id], [host, id]])
            for item in homodb[q]['genes']:
                hostid = str(item[0])
                if hostid in allspp:
                    hostgeneid = str(item[1])
                    # append the original species entrez to entrezdict[0] for
                    # the current species
                    entrezdict[hostid][0].append(q)
                    # append the current species entrez to entrezdict[1]
                    # for the current species
                    entrezdict[hostid][1].append(hostgeneid)
            i += 1

        # sort everything into pandas dataframes
        # for the original species, ensemblg IDs are needed as these are
        # use as the standard in the output database
        # for the other species, gene symbols are needed as these are used
        # by the "mine" databases mined by the "DataMineAnnotation" object.
        # The following code is to translate all the IDs into the appropriate
        # types
        for host in entrezdict:
            # zip original species entrez and other species entrez
            R = list(zip(entrezdict[host][0], entrezdict[host][1]))
            # transalate original species entrez into ensemblg
            R = pd.DataFrame(R, columns=[
                'symbol_%s' % self.ohost, 'entrez_%s' % host])
            R = self.translate(R, 'symbol_%s' % self.ohost, 'ensemblg')
            # download the symbol annotations for the other species
            Sym = SymbolAnnotation(self.datasource, self.outdb,
                                   host=host, shost=names[host])
            Sym.download(entrezdict[host][1], ['symbol'], 'entrezgene',
                         species=host)
            Sym.parse()
            # coerce gene symbols and entrez ids for the other species into a
            # pandas dataframe
            df = Sym.dataset['symbol']

            df = [((key, df[key])) for key in df]
            df = pd.DataFrame(df, columns=['entrez_%s' % host,
                                           'symbol_%s' % host])
            R['entrez_%s' % host] = R['entrez_%s' % host].astype(int)
            df['entrez_%s' % host] = df['entrez_%s' % host].astype(int)
            # merge the two dataframes to generate a translation from original
            # species ensemblg ID to other species gene symbol
            df = df.merge(R,
                          left_on='entrez_%s' % host,
                          right_on='entrez_%s' % host)
            df = df.drop('entrez_%s' % host, 1)

            if not os.path.exists("ensemblg2symbol_%s$geneid" % host):
                self.resultspd['ensemblg2symbol_%s$geneid' % host] = df
                self.storeTranslation(df['ensemblg'],
                                      df['symbol_%s' % host],
                                      ['ensemblg', 'symbol_%s' % host],
                                      'ensemblg2symbol_%s.tsv' % host)


class DataMineAnnotation(APIAnnotation):
    '''
    Downloads data from the standardised "mine" databases which provide
    many annotations for a number of species.
    At least humanmine, mousemine, ratmine, xenmine (xenopus),
    flymine (Drosophila), zebrafishmine and wormmine (C. elegans)
    exist.
    Annotations from these databases can easily be added below.
    The easiest way to find the appropriate views and constraints is using the
    QueryBuilder on the various websites.
    '''

    def __init__(self, prefix, source, outdb, chost, ind, views, constraints,
                 hostid, ohost):
        '''
        views - which fields to download from the database
        constraints - filters on the data to download e.g. list of genes
        these are very specific to the species and the data downloaded so
        should be provided as subclasses below
        host - some databases require a host name in the query, formatted
        e.g. as M. musculus, H. sapiens
        ind - regardless of the order the "views" are listed in, the
        columns in the results will be in a standard order.  ind is the
        column number (0 indexed) which contains the annotation ID.
        '''
        APIAnnotation.__init__(self, prefix, source, outdb, ohost=ohost)
        self.chost = chost
        self.ind = ind
        self.views = views
        self.hostid = hostid
        self.constraints = constraints

    def test(self):
        '''
        Tests the HumanMine API
        Look up symbol for APOBEC3G, should return APOBEC3G.
        '''
        service = SS('http://www.humanmine.org/humanmine/service')
        query = service.new_query("Gene")
        query.add_view("symbol")
        query.add_constraint("Gene", "LOOKUP", "APOBEC3G", code="A")
        for row in query.rows():
            symbol = row['symbol']
        if symbol == "APOBEC3G":
            return 1
        else:
            return 0

    def download(self, genes, fields, scope=None, species=None):
        '''
        Retrives the data depending on self.constraints and self.view
        '''
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
            service = SS(self.datasource)
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
        '''
        Parses the data downloaded using self.download()
        Gene symbols in mouse are translated into ensemblg ids in the original
        species.
        '''
        z = self.dataset
        symbols = []
        terms = []
        other = []
        # ind specifies which column of self.dataset contains the term ID
        # this is used to put the columns into the same order as the column
        # names in the output
        ind = self.ind
        for k in z:
            symbols.append(k[0])
            terms.append(k[ind + 1])
            L = list(k)
            L.remove(k[0])
            L.remove(k[ind + 1])
            k2 = [k[ind + 1]]
            k2 += L
            other.append(tuple(k2))

        # translate symbols into ensemblg IDs in the original species
        dfs = pd.DataFrame(list(zip(symbols, terms)),
                           columns=['symbol_%s' % self.ohost, 'term'])
        dfs = self.translate(dfs, 'symbol_%s' % self.ohost, 'ensemblg')
        # use views to generate column names
        views = self.views
        cols = [v.replace(".", "_") for v in views[0:ind] + views[(ind + 1):]]
        cols = ['term'] + cols
        other = pd.DataFrame(other, columns=cols)
        # store results in the resultspd dictionary to load into the
        # database using self.loadPDTables()
        self.resultspd["%s$details" % self.prefix] = other
        self.resultspd["ensemblg2%s$annot" % self.prefix] = dfs


class MGIAnnotation(DataMineAnnotation):
    '''
    Provides the required options to download, parse and load
    mouse phenotype data using the DataMineAnnotation methods above,
    identified using the mousemine QueryBuilder tool.
    '''

    def __init__(self, source, outdb, ohost):
        prefix = "mgi"
        chost = "Mus_musculus"
        ind = 3
        views = ["ontologyAnnotations.ontologyTerm.description",
                 "ontologyAnnotations.ontologyTerm.name",
                 "ontologyAnnotations.ontologyTerm.namespace",
                 "ontologyAnnotations.ontologyTerm.identifier"]
        constraints = [
            ("Gene.ontologyAnnotations.ontologyTerm.namespace=\
            MPheno.ontology")]
        hostid = 'M. musculus'
        DataMineAnnotation.__init__(self, prefix, source, outdb,
                                    chost, ind, views, constraints, hostid,
                                    ohost=ohost)


class MousePathwayAnnotation(DataMineAnnotation):
    '''
    Provides the required options to load mouse pathway data using
    the DataMineAnnotation methods above, identified using the
    mousemine QueryBuilder tool.
    '''

    def __init__(self, source, outdb, ohost):
        prefix = "mousepathway"
        chost = "Mus_musculus"
        ind = 0
        views = ["pathways.identifier", "pathways.name"]
        constraints = []
        hostid = 'M. musculus'
        DataMineAnnotation.__init__(self, prefix, source, outdb,
                                    chost, ind, views, constraints, hostid,
                                    ohost=ohost)


class HPOAnnotation(DataMineAnnotation):
    '''
    Provides the required options to load human phenotype ontology
    data using the DataMineAnnotation methods above, identifed using
    the humanmine QueryBuilder tool.
    '''

    def __init__(self, source, outdb, ohost):
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
                                    chost, ind, views, constraints, hostid,
                                    ohost=ohost)


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
        '''
        Tests that it is possibly to connect to the URL.
        '''
        try:
            urlopen('http://purl.obolibrary.org/obo/go.owl',
                    timeout=1)
            return True
        except urllib.error.URLError:
            pass

    def download(self, genes=None, fields=None, scope=None, species=None):
        '''
        download an up to date ontology file, parse the xml data into a
        Python "ElementTree" and delete the ontology file.
        '''
        ontologyfile = P.getTempFilename(".")
        os.system("wget -O %s %s" % (ontologyfile, self.datasource))
        tree = ET.parse(ontologyfile)
        os.remove(ontologyfile)
        self.dataset = tree

    def parse(self):
        '''
        Traverse the tree and identify the classes and subclasses
        '''
        # The ns dictionary links "prefixes" of xml tags to the full tag
        # names, see here for more info -
        # https://docs.python.org/2/library/xml.etree.elementtree.html
        ns = {'owl': 'http://www.w3.org/2002/07/owl#',
              'rdfs': 'http://www.w3.org/2000/01/rdf-schema#'}
        tree = self.dataset
        root = tree.getroot()
        r = dict()
        for aclass in root.findall('owl:Class', ns):
            nam = list(aclass.attrib.values())[0].split("/")[-1]
            for subclassof in aclass.findall('rdfs:subClassOf', ns):
                if nam not in r:
                    r[nam] = set()
                if len(list(subclassof.attrib.values())) == 1:
                    r[nam].add(list(subclassof.attrib.values())
                               [0].split("/")[-1])

        klist = []
        vlist = []
        for item in list(r.items()):
            key = item[0]
            vals = item[1]
            for v in vals:
                klist.append(key.replace("_", ":"))
                vlist.append(v.replace("_", ":"))
        self.resultsz['%s$ont' % self.prefix] = [list(zip(klist, vlist)),
                                                 ['id', 'subclass_of'],
                                                 ['varchar(250)',
                                                  'varchar(250)']]


@cluster_runnable
def runall(obj, genelist=None, fields=None, scope='symbol', species=None):
    '''
    Runs object.runall on the cluster
    '''
    obj.runall(genelist, fields, scope, species)


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
    tables = [tab[0] for tab in tables]
    D = dict()
    for tab in tables:
        ttype = tab.split("$")[1]
        D.setdefault(ttype, [])
        D[ttype].append(tab)
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
def MakeSubDBs(infile, newdb, idtype, db):
    '''
    Subsets the database for user specified lists of genes - generates a new
    database for each input file containing annotations for the genes in the
    input only - allows the user to quickly see details of their genes.
    '''
    dbh = sqlite3.connect(newdb)
    genelist = readGeneList(infile)
    if idtype == "ensemblg":
        tgenelist = genelist
    else:
        translation = "ensemblg2%s.tsv" % idtype
        trans = pd.read_csv(translation, sep="\t")
        translation = trans['ensemblg'][trans[idtype].isin(genelist)]
        tgenelist = list(translation)
    tabs = getTables(db)

    for tab in tabs['geneid'] + tabs['other']:
        if 'ensemblg' in tab:
            res = readDBTable(db, tab)
            cnames = getDBColumnNames(db, tab)
            res = pd.DataFrame(res, columns=cnames)
            subtab = res[res['ensemblg'].isin(tgenelist)]
            tname = tab.split("$")[0].replace("ensemblg2", "")
            if idtype != 'ensemblg':
                subtab = subtab.merge(trans, 'left', left_on='ensemblg',
                                      right_on='ensemblg')
                subtab = subtab.drop('ensemblg', 1)
            subtab.to_sql(tname, dbh, index=False, if_exists='replace')

    for tab in tabs['annot']:
        if 'ensemblg' in tab:
            annot = readDBTable(db, tab)
            cnames = getDBColumnNames(db, tab)
            annot = pd.DataFrame(annot, columns=cnames)
            detailstab = tab.replace(
                "annot", "details").replace("ensemblg2", "")
            details = readDBTable(db, detailstab)
            cnames = getDBColumnNames(db, detailstab)
            details = pd.DataFrame(details, columns=cnames)

            mergeon = set(annot.columns)
            mergeon.remove('ensemblg')
            mergeon = list(mergeon)[0]
            subannot = annot[annot['ensemblg'].isin(tgenelist)]
            details = details[details[mergeon].isin(subannot[mergeon])]
            tname = tab.split("$")[0].replace("ensemblg2", "")
            if idtype != 'ensemblg':
                subannot = subannot.merge(
                    trans, 'left', left_on='ensemblg', right_on='ensemblg')
                subannot = subannot.drop('ensemblg', 1)
            subannot.to_sql("%s_annotations" % (tname), dbh,
                            index=False, if_exists='replace')
            details.to_sql("%s_details" % (tname), dbh,
                           index=False, if_exists='replace')
