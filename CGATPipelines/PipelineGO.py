"""PipelineGO - Tasks for a GO analysis
============================================

Reference
---------


"""
import re
import glob
import os
import collections
import sqlite3

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.Stats as Stats
import CGAT.IOTools as IOTools
import CGAT.CSV as CSV

# set from calling module
PARAMS = {}


def createGOFromENSEMBL(infile, outfile):
    """get GO assignments from ENSEMBL

    Download GO assignments from the ENSEMBL database and store in
    file.

    Configuration
    -------------
    go_host
       ENSEMBL mysql host name, e.g., ensembldb.ensembl.org
    go_database
       ENSEMBL database, e.g., homo_sapiens_core_75_37
    go_port
      ENSEMBL port to use, e.g., 5306

    Arguments
    ---------
    infile : string
        Unused
    outfile : string
        Output filename

    """

    job_memory = "5G"
    statement = '''
    cgat runGO
    --filename-dump=%(outfile)s
    --database-host=%(go_host)s
    --database-user=anonymous
    --database-name=%(go_database)s
    --database-port=%(go_port)i
    > %(outfile)s.log
    '''

    P.run()


def createGOFromGeneOntology(infile, outfile):
    """get GO assignments from Geneontology.org

    GO terms are mapped to ensembl gene names via uniprot identifiers.

    Configuration
    -------------
    geneontology_file
       Filename on geneontology database, e.g.,
       gene_association.goa_human.gz
    database_name
       Pipeline database name

    Arguments
    ---------
    infile : string
        Unused
    outfile : string
        Output filename
    """

    filename = os.path.join(os.path.dirname(outfile), "geneontology.goa.gz")
    if not os.path.exists(filename):
        statement = '''
    wget -O %(filename)s http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/%(go_geneontology_file)s?rev=HEAD
    '''

        P.run()

    # see http://www.geneontology.org/gene-associations/readme/goa.README
    Data = collections.namedtuple(
        "Data",
        "db db_object_id db_object_symbol qualifier goid dbreference evidence "
        " with_id aspect "
        " db_object_name synonym db_object_type "
        " taxon_id date assigned_by "
        " annotation_extension"
        " gene_product_form_id")

    dbh = sqlite3.connect(PARAMS["database_name"])
    cc = dbh.cursor()
    map_uniprot2ensembl = dict(
        cc.execute("SELECT DISTINCT gene_name, gene_id FROM transcript_info").fetchall())
    map_goid2description = dict(
        cc.execute("SELECT DISTINCT go_id, description FROM go_assignments").fetchall())

    aspect2name = {"P": "biol_process",
                   "F": "mol_function",
                   "C": "cell_location"}

    c = E.Counter()
    found_uniprot, found_genes, notfound_uniprot = set(), set(), set()
    outf = IOTools.openFile(outfile, "w")
    outf.write("go_type\tgene_id\tgo_id\tdescription\tevidence\n")
    for line in IOTools.openFile(filename):
        if line.startswith("!"):
            continue
        c.input += 1
        data = Data._make(line[:-1].split("\t"))

        if data.db_object_symbol in map_uniprot2ensembl:
            gene_id = map_uniprot2ensembl[data.db_object_symbol]
            found_uniprot.add(data.db_object_symbol)
            found_genes.add(gene_id)
            outf.write("%s\t%s\t%s\t%s\t%s\n" %
                       (aspect2name[data.aspect],
                        gene_id,
                        data.goid,
                        map_goid2description.get(data.goid, ""),
                        data.evidence))
            c.output += 1

        else:
            c.notfound += 1
            notfound_uniprot.add(data.db_object_symbol)

    c.found_genes = len(found_genes)
    c.found_uniprot = len(found_uniprot)
    c.notfound_uniprot = len(notfound_uniprot)

    E.info("%s" % str(c))
    E.info("not found=%s" % str(notfound_uniprot))
    outf.close()

############################################################
############################################################
############################################################
# get GO descriptions
############################################################


def imputeGO(infile_go, infile_paths, outfile):
    """impute GO accessions.

    Output a list of gene-to-GO associations for genes that includes
    ancestral terms.

    Arguments
    ---------
    infile_go : string
        Filename with gene-to-GO assocations for genes
    infile_paths : string
        Filename with paths of term to ancestor (see go2fmt.pl).
    outfile : string
         Output filename

    """

    c = E.Counter()

    term2ancestors = collections.defaultdict(set)
    with IOTools.openFile(infile_paths) as inf:
        for line in inf:
            parts = line[:-1].split()
            term = parts[0]
            ancestors = [parts[x] for x in range(2, len(parts), 2)]
            # there can be multiple paths
            term2ancestors[term].update(ancestors)

    goid2description = {}
    gene2goids = collections.defaultdict(list)
    goid2type = {}
    with IOTools.openFile(infile_go) as inf:
        for line in inf:
            if line.startswith("go_type"):
                continue
            go_type, gene_id, goid, description, evidence = line[
                :-1].split("\t")
            gene2goids[gene_id].append(goid)
            goid2description[goid] = description
            goid2type[goid] = go_type

    outf = IOTools.openFile(outfile, "w ")
    for gene_id, in_goids in gene2goids.items():
        c.genes += 1
        out_goids = set(in_goids)
        for goid in in_goids:
            out_goids.update(term2ancestors[goid])
        if len(in_goids) != len(out_goids):
            c.increased += 1
        else:
            c.complete += 1

        for goid in out_goids:
            outf.write("\t".join(
                (goid2type.get(goid, ""), gene_id, goid,
                 goid2description.get(goid, ""), "NA")) + "\n")
            c.assocations += 1

    outf.close()

    E.info("%s" % str(c))


def buildGOPaths(infile, outfile):
    """create paths from derived to ancestral terms.

    The output of this function maps all terms inside an ontology to
    all its ancestors.

    Arguments
    ---------
    infile : string
       An ontology (.obo) file
    outfile : string
       Output file with paths of terms to root.

    """

    statement = '''
    go2fmt.pl -w pathlist %(infile)s > %(outfile)s
    '''
    P.run()


def buildGOTable(infile, outfile):
    """convert ontology to tab-separated table.

    Arguments
    ---------
    infile : string
        An ontology (.obo) file
    outfile : string
        Output file mapping GO terms to their description,
        long description and text.
    """

    statement = '''
    echo -e "go_id\\tdescription\\tlong_description\\ttext\\n" > %(outfile)s;
    go2fmt.pl -w tbl %(infile)s >> %(outfile)s
    '''
    P.run()


def getGODescriptions(infile):
    '''build dictionary mapping GOids to types and descriptions.

    Arguments
    ---------
    infile : string
        Filename of table with GO assignments

    Returns
    -------
    mapping : dict
        Dictionary mapping GOid to GOtype and GOdescription.
    '''

    with IOTools.openFile(infile) as inf:
        fields, table = CSV.readTable(inf, as_rows=False)

    return dict([(y, (x, z)) for x, y, z in zip(
        table[fields.index("go_type")],
        table[fields.index("go_id")],
        table[fields.index("description")])])


def createGOSlimFromENSEMBL(infile, outfile):
    """build GO SLIM assignments.

    This method downloads a GOSlim specification
    from ``go_url_goslim`` and applies it to
    a table of GO assignments.

    Config
    ------
    go_url_goslim
       GOSlim definition, e.g.,
       http://www.geneontology.org/ontology/subsets/goslim_generic.obo
    go_url_ontology
       GO Ontology location, e.g.,
       http://www.geneontology.org/ontology/gene_ontology.obo

    Arguments
    ---------
    infile : string
        Filename with GO assignments
    outfile : strinng
        Output filename with GOSlim assignments.

    """

    dirname = os.path.dirname(outfile)

    E.info("downloading GOSlim specification from %s" %
           PARAMS["go_url_goslim"])
    goslim_fn = os.path.join(dirname, "goslim.obo")
    statement = '''wget %(go_url_goslim)s
    --output-document=%(goslim_fn)s'''
    P.run()

    E.info("downloading GO ontology from %s" %
           PARAMS["go_url_ontology"])
    ontology_fn = os.path.join(dirname, "go_onotology.obo")
    statement = '''wget %(go_url_ontology)s
    --output-document=%(ontology_fn)s'''
    P.run()

    E.info("mapping GO to GOSlim")
    job_memory = "5G"
    statement = '''
    map2slim -outmap %(outfile)s.map
    %(goslim_fn)s
    %(ontology_fn)s
    '''
    P.run()

    job_memory = "5G"
    statement = '''
    zcat < %(infile)s
    | cgat runGO
    --go2goslim
    --filename-ontology=%(ontology_fn)s
    --slims=%(outfile)s.map
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


def runGOFromFiles(outfile,
                   outdir,
                   fg_file,
                   bg_file=None,
                   go_file=None,
                   ontology_file=None,
                   samples=None,
                   minimum_counts=0,
                   pairs=False,
                   gene2name=None):
    """check for GO enrichment.

    The gene lists are supplied by files.
    This method is a wrapper for `runGO.py`.

    Arguments
    ---------
    outfile : string
        Output filename
    outdir : string
        Output directory for auxiliary files
    fg_file : string
        Gene list of foreground.
    bg_file : string
        Gene list for background. If None, all genes
        with GO annotations are used as background.
    go_file : string
        Filename with Gene-to-GO assignments
    ontology_file : string
        Filename with ontology information.
    samples : int
        Number of samples for empirical FDR. If not given, use
        BH FDR.
    minimum_counts : int
        Minimum number of observations in a GO category
        required in order for using it.
    pairs : bool
       If True, each category for each pair of gene sets will
       be tested for differential enrichment.
    gene2name : string
        Filename with a gene-to-genename information.
    """

    if ontology_file is None:
        ontology_file = PARAMS.get("go_ontology", None)

    options = []
    if ontology_file:
        options.append("--filename-ontology=%(ontology_file)s" % locals())

    if bg_file is not None:
        options.append("--background-tsv-file=%(bg_file)s" % locals())

    if samples is not None:
        options.append("--fdr")
        options.append("--sample-size=%(samples)i" % locals())
        options.append("--fdr-method=empirical")
    else:
        options.append("--fdr")
        options.append("--fdr-method=BH")

    if pairs:
        options.append("--pairwise")

    if gene2name:
        options.append("--gene2name-map-tsv-file=%s" % gene2name)

    options = " ".join(options)
    statement = '''
    cgat runGO 
    --filename-input=%(go_file)s
    --genes-tsv-file=%(fg_file)s
    --output-filename-pattern='%(outdir)s/%%(set)s.%%(go)s.%%(section)s'
    --min-counts=%(minimum_counts)i
    --log=%(outfile)s.log
    %(options)s
    > %(outfile)s'''

    P.run()


def runGOFromDatabase(outfile,
                      outdir,
                      statement_fg,
                      statement_bg,
                      go_file,
                      ontology_file=None,
                      samples=1000):
    """check for GO enrichment.

    Gene lists are extracted from a database.
    This method is a wrapper for `runGO.py`.

    Arguments
    ---------
    outfile : string
        Output filename
    outdir : string
        Output directory for auxiliary files
    statement_fg : string
        SQL statement to select genes of foreground set.
    statement_bg : string
        SQL statement to select genes in background set.
    go_file : string
        Filename with Gene-to-GO assignments
    ontology_file : string
        Filename with ontology information.
    samples : int
        Number of samples for empirical FDR. If not given, use
        BH FDR.
    """

    dbhandle = sqlite3.connect(PARAMS["database_name"])

    cc = dbhandle.cursor()
    fg = set([x[0] for x in cc.execute(statement_fg).fetchall()])
    bg = set([x[0] for x in cc.execute(statement_bg).fetchall()])

    if len(fg) == 0:
        P.touch(outfile)
        return

    fg_file = os.path.join(outdir, "foreground")
    bg_file = os.path.join(outdir, "background")
    outf = open(fg_file, "w")
    outf.write("\n".join(map(str, fg)) + "\n")
    outf.close()
    outf = open(bg_file, "w")
    outf.write("\n".join(map(str, bg)) + "\n")
    outf.close()

    runGOFromFiles(outfile, outdir,
                   fg_file, bg_file,
                   go_file,
                   ontology_file=ontology_file,
                   samples=samples)


def loadGO(infile, outfile, tablename):
    """import GO results into individual tables.

    This method concatenates all the results from
    a GO analysis and uploads into a single table.

    """

    indir = infile + ".dir"

    if not os.path.exists(indir):
        P.touch(outfile)
        return

    load_statement = P.build_load_statement(
        tablename=tablename,
        options="--allow-empty-file "
        "--add-index=category "
        "--add-index=goid ")

    statement = '''
    python %(toolsdir)s/cat_tables.py %(indir)s/*.overall
    | %(load_statement)s
    > %(outfile)s
    '''
    P.run()


def loadGOs(infiles, outfile, tablename):
    '''import GO results into a single table.

    This method also computes a global QValue over all
    tracks, genesets and annotation sets.

    Arguments
    ---------
    infiles : string
       Output files of several runGO analyses
    outfile : string
       Output filename, contains log information
    tablename : string
       Table name for storing results.
    '''

    header = False

    tempf1 = P.getTempFile()

    pvalues = []

    for infile in infiles:
        indir = infile + ".dir"

        if not os.path.exists(indir):
            continue

        track, geneset, annotationset = re.search(
            "^(\S+)_vs_(\S+)\.(\S+)", infile).groups()

        for filename in glob.glob(os.path.join(indir, "*.overall")):
            for line in open(filename, "r"):
                if line.startswith("#"):
                    continue
                data = line[:-1].split("\t")
                if line.startswith("code"):
                    if header:
                        continue
                    tempf1.write("track\tgeneset\tannotationset\t%s" % line)
                    header = True
                    assert data[10] == "pover" and data[
                        11] == "punder", "format error, expected pover-punder, got %s-%s" % (data[10], data[11])
                    continue
                tempf1.write("%s\t%s\t%s\t%s" %
                             (track, geneset, annotationset, line))
                pvalues.append(min(float(data[10]), float(data[11])))

    tempf1.close()

    E.info("analysing %i pvalues" % len(pvalues))
    fdr = Stats.doFDR(pvalues)
    E.info("got %i qvalues" % len(fdr.mQValues))
    qvalues = ["global_qvalue"] + fdr.mQValues

    tempf2 = P.getTempFile()

    for line, qvalue in zip(open(tempf1.name, "r"), qvalues):
        tempf2.write("%s\t%s\n" % (line[:-1], str(qvalue)))

    tempf2.close()

    P.load(tempf2.name, outfile,
           tablename=tablename,
           options="--allow-empty-file "
           "--add-index=category "
           "--add-index=track,geneset,annotationset "
           "--add-index=geneset "
           "--add-index=annotationset "
           "--add-index=goid ")

    os.unlink(tempf1.name)
    os.unlink(tempf2.name)
