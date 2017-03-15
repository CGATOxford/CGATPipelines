"""
====================
Gene Info Pipeline
====================

:Author: Katy Brown
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
========

This pipeline generates a database of annotations mapped to ensembl gene ids
for all current genes in Entrez Gene for a chosen species.
The annotations are downloaded from various APIs.
This database provides annotation information about the genes and is used
as an input for pipeline_enrichment.py

annotate - load all annotations specified in pipeline.ini
full - load all annotations and if specified in pipeline.ini, generate
a database for each list of genes in genelists.dir containing annotations
for genes in this list

Usage and Inputs
================
The pipeline requires a configured :file:`pipeline.ini` file.
No other inputs are required.


Output
======
The major output is a database, named as specified in the pipeline.ini.
All annotations are mapped to ensemblg ids.

Tables in the database are named as:
- xxxx$geneid - translations from ensemblg to another type of gene ID
- xxxx$annot - ensemblg IDs and a corresponding annotation ID.  These
  annotations can be used later in enrichment analysis
- xxxx$details - annotation IDs and any columns containing metadata about
  this annotation
- xxxx$ont - terms from a hierarchical ontology and the direct parents of that
  term
- xxxx$other - miscellanous annotations not to be used for enrichment
  analysis


Annotation
==========
Annotations are generated using APIAnnotation objects in the module
file - these provide the following methods.
test() - connectivity test for the API
download() - downloads data from the API
parse() - parses the data into a useable format - either a list
of tuples or a pandas dataframe
loadPDTables - loads data stored in a dictionary of pandas dataframes
loadZippedTables - loads data stored in a dictionary of lists of tuples

Annotations are currently generated using four APIs and different subclasses
of APIAnnotation are used accordingly.
Subclasses marked with (x) are always generated, all others are optional.

EntrezAnnotation
Downloaded from the Entrez API via the Entrez BioPython module.

- (x) EntrezGeneAnnotation - used to generate the initial gene list only
- EntrezTaxonomyAnnotation - matches species to taxonomic info - used
  to retrieve scientific names of species

MyGeneInfoAnnotation
mygene.info provides many regularly updated annotations via its API

- (x) SymbolAnnotation - translates Ensembl gene IDs to gene symbols
- (x) EnsemblAnnotation - translates Entrez IDs to Ensembl Gene, Ensembl
   Transcript and Ensembl Protein IDs
- GoAnnotation - gene ontology annotation (any or all of BP, MF, CC)
- PathwayAnnotation - pathway annotation (any or all of kegg, humancyc,
  biocarta, mousecyc, netpath, pharmgkb, pid, reactome, smpdb, wikipathways,
  yeastcyc)
- HomologeneAnnotation - gene symbols for homologous genes in other species
  as listed in the pipeline.ini

DataMineAnnotation
Communicates with the "mine" databases - humanmine.org, mousemine.org,
ratmine.org etc. which provide many annotations for these species

- HPOAnnotation - human phenotype ontology annotation
- MGIAnnotation - mouse phenotype annotation (based on homologous genes)

OntologyAnnotation
Communicates with the obolibrary (obo foundry) API for hierarchical
ontology annotations.
Can be used to download and parse any OWL formatted ontology available
on this site.

"""

from ruffus import *
import CGATPipelines.Pipeline as P
import os
import sys
import CGATPipelines.PipelineGeneInfo as PipelineGeneInfo
import CGAT.IOTools as IOTools
import pandas as pd


PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# pick a pathway to use to check pathway annotation has run
example_pw = PARAMS['my_gene_info_pathway'].split(",")[0]
if example_pw == "all":
    example_pw = 'kegg'
example_homolo = str(PARAMS['my_gene_info_homologene']).split(",")
if len(example_homolo) == 0 or example_homolo[0] == 'all':
    example_homolo = 10090
else:
    example_homolo = example_homolo[0]

# get the list of annotations to be downloaded from my gene info
mgiannotations = PARAMS['my_gene_info_annotations'].split(",")


@originate('allgenes.tsv')
def GetAndTranslateAllGenes(outfile):
    '''
    This step is required.
    1. All Entrez gene IDs are downloaded from entrez gene.
    2. Corresponding ensembl gene, ensembl transcript and ensembl protein
       IDs are downloaded from mygene.info
    3. Corresponding gene symbols are downloaded from mygene.info
    4. These are loaded into the database
    5. A list of all gene Entrez IDs is stored as 'allgenes.tsv

    Tables:
    ensemblg2entrez$geneid - ensemblg to entrez ID
    ensemblg2ensemblt$other - ensemblg to ensembl transcript
    ensemblg2ensemblp$other - ensemblg to ensembl protein
    ensemblg2symbol_xxx$geneid - ensemblg to symbol in species xxx
    '''
    GeneAnnot = PipelineGeneInfo.EntrezGeneAnnotation(
        PARAMS['db_name'], PARAMS['entrez_email'])
    if PARAMS['test'] == 1:
        entrezgenelist = GeneAnnot.download_all(PARAMS['entrez_host'],
                                                count=100)
    else:
        entrezgenelist = GeneAnnot.download_all(PARAMS['entrez_host'])

    # Generate a SymbolAnnotation object
    Sym = PipelineGeneInfo.SymbolAnnotation(PARAMS['my_gene_info_source'],
                                            PARAMS['db_name'],
                                            PARAMS['entrez_host'],
                                            PARAMS['entrez_sciname'])

    # Get Symbol Annotations
    PipelineGeneInfo.runall(Sym, entrezgenelist, ['symbol'],
                            scope='entrezgene', species=PARAMS['entrez_host'],
                            submit=True)

    genesymbols = list(pd.read_csv("entrez2symbol_%s.tsv" % PARAMS[
        'entrez_host'], sep="\t")['symbol_%s' % PARAMS['entrez_host']])

    # Generate an EnsemblAnnotation object
    Ens = PipelineGeneInfo.EnsemblAnnotation(PARAMS['my_gene_info_source'],
                                             PARAMS['db_name'],
                                             PARAMS['entrez_host'])
    # Get Ensembl annotations
    PipelineGeneInfo.runall(Ens, genesymbols, ['ensembl'], scope="symbol",
                            species=PARAMS['entrez_host'], submit=True)

    # Make output gene list
    outf = IOTools.openFile(outfile, "w")
    for gene in genesymbols:
        outf.write("%s\n" % gene)
    outf.close()


@active_if('go' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2go\$annot.load')
def AnnotateWithGO(infile, outfile):
    '''
    Annotates all genes in allgenes.tsv with GO ontology terms using
    information from mygene.info
    Tables:
    ensemblg2go$annot- ensemblg to go ID
    go$details - go ID to details of go term
    go$ont - go ID to parent go IDs
    '''
    genelist = PipelineGeneInfo.readGeneList(infile)
    # Generate a GoAnnotation object with details from mygene.info
    GO = PipelineGeneInfo.GoAnnotation(PARAMS['my_gene_info_source'],
                                       PARAMS['db_name'],
                                       PARAMS['my_gene_info_go'],
                                       PARAMS['entrez_host'])
    PipelineGeneInfo.runall(GO, genelist, ['go'],
                            species=PARAMS['entrez_host'], submit=True)

    # Get the GO hierarcical ontology from OBO foundry
    ont = PipelineGeneInfo.OntologyAnnotation('go',
                                              PARAMS['my_gene_info_goont'],
                                              PARAMS['db_name'])
    PipelineGeneInfo.runall(ont, genelist, species=PARAMS['entrez_host'],
                            submit=True)


@active_if('pathway' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2%s\$annot.load' % example_pw)
def AnnotateWithPathway(infile, outfile):
    '''
    Annotates all genes in allgenes.tsv with pathway details, either
    for all pathway databases available via mygene.info or those
    specified in the pipeline.ini
    Tables:
    ensemblg2xxx$annot - ensemblg to ID in pathway database
    xxx$details - pathway database ID to pathway details
    '''
    genelist = PipelineGeneInfo.readGeneList(infile)
    PW = PipelineGeneInfo.PathwayAnnotation(PARAMS['my_gene_info_source'],
                                            PARAMS['db_name'],
                                            PARAMS['my_gene_info_pathway'],
                                            PARAMS['entrez_host'])
    PipelineGeneInfo.runall(PW, genelist,
                            ['pathway'], species=PARAMS['entrez_host'],
                            submit=True)


@active_if('homologene' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2symbol_%s\$annot.load' % example_homolo)
def AnnotateWithHomologene(infile, outfile):
    '''
    Annotates all genes in allgenes.tsv with homologous gene symbols from
    either a list of species provided in the pipeline.ini or all species
    available in homologene via mygene.info
    Tables:
    ensemblg2symbol_xxx$geneid - ensemblg in original species to symbol in xxx
    '''
    genelist = PipelineGeneInfo.readGeneList(infile)
    HG = PipelineGeneInfo.HomologeneAnnotation(PARAMS['my_gene_info_source'],
                                               PARAMS['db_name'],
                                               PARAMS[
                                                   'my_gene_info_homologene'],
                                               PARAMS['entrez_host'],
                                               PARAMS['entrez_email'])
    PipelineGeneInfo.runall(HG, genelist, ['homologene'],
                            species=PARAMS['entrez_host'], submit=True)


@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_mousepathway']) == 1)
@transform('ensemblg2symbol_10090$geneid.load', suffix(".load"),
           'ensemblg2mousepathway\$annot.load')
def AnnotateWithMousePathway(infile, outfile):
    '''
    Uses the list of mouse gene symbols generated using homologene above
    to annotate mouse pathways provided at mousemine.org
    Tables:
    ensemblg2mousepathway$annot - original host ensemblg to
                                  mouse pathway ID
    mousepathway$details - mouse pathway ID to mouse pathway details
    '''
    genelist = list(set(PipelineGeneInfo.getSymbols(infile)))
    MP = PipelineGeneInfo.MousePathwayAnnotation(
        PARAMS['homologues_mousemine'],
        PARAMS['db_name'], ohost=PARAMS['entrez_host'])
    PipelineGeneInfo.runall(MP, genelist, submit=True)


@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_mgi']) == 1)
@transform('ensemblg2symbol_10090$geneid.load', suffix(".load"),
           'ensemblg2mgi\$annot.load')
def AnnotateWithMGI(infile, outfile):
    '''
    Uses the list of mouse gene symbols generated using homologene above
    to annotate mouse phenotypes provided through MGI at mousemine.org
    Tables:
    ensemblg2mgi$annot - original host ensemblg to mouse phenotype ID
    mgi$details - mouse phenotype ID to mouse phenotype details
    '''
    genelist = list(set(PipelineGeneInfo.getSymbols(infile)))
    MGI = PipelineGeneInfo.MGIAnnotation(
        PARAMS['homologues_mousemine'],
        PARAMS['db_name'], ohost=PARAMS['entrez_host'])
    PipelineGeneInfo.runall(MGI, genelist, submit=True)


@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_hpo']) == 1)
@transform('ensemblg2symbol_9606$geneid.load', suffix(".load"),
           'ensemblg2hpo\$annot.load')
def AnnotateWithHPO(infile, outfile):
    '''
    Uses the list of human gene symbols generated using homologene above
    to annotate human phenotypes provided through HPO at humanmine.org
    Tables:
    ensemblg2hpo$annot - original host ensemblg to human phenotype ID
    hpo$details - human phenotype ID to human phenotype details
    '''
    genelist = list(set(PipelineGeneInfo.getSymbols(infile)))
    HPO = PipelineGeneInfo.HPOAnnotation(
        PARAMS['homologues_humanmine'],
        PARAMS['db_name'], PARAMS['entrez_host'])
    PipelineGeneInfo.runall(HPO, genelist, submit=True)
    ont = PipelineGeneInfo.OntologyAnnotation('hpo',
                                              PARAMS['homologues_hpoont'],
                                              PARAMS['db_name'])
    PipelineGeneInfo.runall(ont, genelist, species=PARAMS['entrez_host'],
                            submit=True)


@follows(AnnotateWithGO)
@follows(AnnotateWithPathway)
@follows(AnnotateWithHomologene)
@follows(AnnotateWithMGI)
@follows(AnnotateWithMousePathway)
@follows(AnnotateWithHPO)
@transform("allgenes.tsv", suffix(".tsv"), '.tsv')
def annotate(infile, outfile):
    pass


@follows(mkdir('genesetdbs.dir'))
@follows(annotate)
@active_if(int(PARAMS['db_subset']) == 1)
@transform("genelists.dir/*.tsv", regex("genelists.dir/(.*).tsv"),
           r"genesetdbs.dir/\1")
def MakeSubDBs(infile, outfile):
    '''
    Takes any lists of genes provided in genesets.dir and makes a database
    in genesetdbs.dir containing only annotations for genes in the list.
    These will have the same gene ID type as the input lists
    and allow the user to quickly see the annotations for their genes
    of interest.
    '''
    PipelineGeneInfo.MakeSubDBs(infile, outfile, PARAMS['db_subsettype'],
                                PARAMS['db_name'], submit=True)


@follows(MakeSubDBs)
def full():
    pass


@follows(MakeSubDBs)
def build_report():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
