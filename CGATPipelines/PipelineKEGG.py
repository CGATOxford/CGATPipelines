'''
PipelineKegg.py - Tasks related to KEGG analysis
================================================

Reference
---------
'''

import CGAT.Experiment as E
from rpy2.robjects import r as R
import CGAT.IOTools as IOTools
import CGAT.Biomart as Biomart
import re


def importKEGGAssignments(outfile, mart, host, biomart_dataset):
    '''import the KEGG annotations from the R KEGG.db annotations
    package.

    .. note::

        Since KEGG is no longer publically available, this is not
        up-to-date and maybe removed from bioconductor in future
        releases

    The table written to outfile has the following columns:
    ``ontology``, ``gene_id``, ``kegg_ID``, ``kegg_name``,
    ``evidence``.

    Arguments
    ---------
    outfile : string
        Output filename in :term:`tsv` format.
    mart : string
        Name of the biomart
    host : string
        Host name of the biomart server
    biomart_dataset : string
        Biomart dataset

    '''

    if not re.match("rnorvegicus|scerevisiae|hsapiens|mmusculus",
                    biomart_dataset):
        E.warn("KEGG.db doesn't map Entrez ids for %s, %s will"
               " likely be empty" % (biomart_dataset, outfile))

    R.library("KEGG.db")

    E.info("getting entrez to ensembl mapping ...")
    entrez2ensembl = Biomart.biomart_iterator(
        ("ensembl_gene_id", "entrezgene"),
        biomart=mart,
        dataset=biomart_dataset,
        host=host,
        path="/biomart/martservice")

    entrez2ensembl = dict((x['entrezgene'],
                           x['ensembl_gene_id'])
                          for x in entrez2ensembl)

    E.info("Done")

    E.info("getting entrez to kegg mapping ... ")
    entrez2path = R('as.list(KEGGEXTID2PATHID)')
    E.info("Done")

    E.info("Getting KEGG names")
    pathnames = R('as.list(KEGGPATHID2NAME)')
    pathid2name = dict(list(zip(pathnames.names, R.unlist(pathnames))))
    E.info("Done")

    outf = IOTools.openFile(outfile, "w")
    outf.write("ontology\tgene_id\tkegg_ID\tkegg_name\tevidence\n")

    # rx2 did not work in rpy2 2.4.2 - workaround uses
    # absolute indices
    for gene_column, gene in enumerate(entrez2path.names):

        try:
            gene = int(gene)
        except ValueError:
            continue

        if gene in entrez2ensembl:
            ensid = entrez2ensembl[gene]

        else:
            continue

        for pathway in entrez2path[gene_column]:
            pathid = re.match("[a-z]+([0-9]+)", pathway).groups()[0]
            pathname = pathid2name[pathid]
            outf.write(
                "\t".join(["kegg", ensid, str(pathway),
                           pathname, "NA"]) + "\n")
