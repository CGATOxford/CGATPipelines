##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline Enrichment
===========================

:Author: Katy Brown
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline calculates if any of a list of "terms" (any list of identifiers
which has been mapped to correspond to a list of genes) are significantly
more often mapped to a gene in the "foreground" - a list of genes identified
in an experiment, e.g. differentially expressed genes, than would be expected
given the "background" - a list of all genes which could have been identified
in the experiment.
For example, if looking for enrichment in the list of genes which are
differentially expressed in diseased vs control kidney samples, the
foreground would be the list of significantly differentially expressed genes
and the background could be, for example, a list of all genes known to be
expressed in the kidney.

# IMPORTANT #
Gene IDs specified more than once in the input gene lists will be filtered
out to leave only unique gene names - there is no weighting according to
how many times the gene is in the input list.


Initially, the working directory should contain:
a directory, "foregrounds.dir"
This should contain a file for each list of foreground genes, with the
suffix .tsv, with one gene per line.  Any gene identifier can be used, but
there needs to be a table in the get_gene_info.py database named
ensemblg2[identifier_type]$annot mapping ensemblg identifiers to this
identifier type.  All foregrounds need to use the same identifier type.
(Currently ensemblg, gene symbol and entrez id are implemented)

optionally, a directory "backgrounds.dir"
This should contain any user specific backgrounds to compare to the
foregrounds.  These should be formatted as for the foregrounds and
should all have the same identifier type to each other (but not necessarily
to the foregrounds).
All foregrounds will automatically be compared to a background of all genes
with annotations of the specified type, so these backgrounds do not need
to be provided.

Gene Identifiers
# IMPORTANT #
There is not a 1:1 relationshup between different types of gene identifier
(i.e. ensembl gene, entrez and hgnc) so you need to decide which
list to use and stick to it.

In the pipeline.ini you need to specify the type of gene ID you want to use.
Foreground and background input lists should use this type of ID.
For the annotation steps of the pipeline ensembl gene IDs are used,
but before the enrichment test the gene IDs are translated back to the
initial type and duplicates are removed.
The "all genes" background will be all identifiers of the specified type
- so if you choose ensemblg it will be all ensembl gene IDs including those
which correspond to the same hgnc ID.

If the input list has the same gene ID twice it will be filtered out,
but if the input list has two gene IDs which correspond to the same gene
it will not be filtered out.


Annotations
Annotations are stored in AnnotationSet objects.
These contain all the annotations associated with a particular
annotation source.

Annotations from a Database
AnnotationSets are predominantly generated from a database using an
AnnotationParser method.
The Database is generated using the pipeline pipeline_geneinfo.py.
This database is required to run pipeline_enrichment

The tables in this database corresponding to annotations are in sets of two
or three, named prefix$annot, prefix$details and prefix$ont, where prefix
is the name of the annotation source.

$annot contains each gene and each term it is annotated to.  There should be
one gene and term per row, genes mapped to multiple terms are
repeated multiple times e.g.
    gene term
    A    1
    A    2
    B    2

$details has one row for each term, the remaining columns are metadata
associated with the term in the annotation source.

$ont contains each term mapped to its parent terms (in an ontology) - if
the annotation source is not hierarchical this is not needed.


Annotations from a Flat File
The user can specify additional annotations to use which are not in the
database as long as a file exists containing all the required information.
Processing of this file requires users to specify a set of
"annotations2annotations.py" options in the pipeline.ini.


AnnotationSets
Annotations are converted to AnnotationSet objects.

An AnnotationSet consists of four dictionaries:
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


Enrichment Testing
Every foreground will be tested for enrichment against every background using
every available annotation.

There are various ways to test for enrichment, these are contained in
EnrichmentTester objects.  The standard is to test for enrichment of each
term individually, however alternatives can be added as subclasses of
EnrichmentTester and specified in the pipeline.ini.

Statistical tests for enrichment are called by the EnrichmentTester class
and specified as StatsTest subclasses.  The test type, multiple testing
correction and signficance threshold are specified in the pipeline.ini.
(Only FisherExactTest currently implemented)


"""
from ruffus import *
from ruffus.combinatorics import *
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import sys
import CGATPipelines.PipelineGSEnrichment as PipelineEnrichment
import CGAT.IOTools as IOTools
import pandas as pd


# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

dbname = PARAMS['db_name']
unmapped = PipelineEnrichment.getUnmapped(PARAMS)
outfilesuffixes = ["_genestoterms.tsv",
                   "_termstogenes.tsv",
                   "_termstodetails.tsv",
                   "_termstoont.tsv"]
unmappedouts = [["annotations.dir/%s%s" % (u, s)
                 for s in outfilesuffixes]
                for u in unmapped]

hpatissues = PARAMS['hpa_tissue'].split(",")
hpatissues = ['clean_backgrounds.dir/%s_hpa_background.tsv' % tissue
              for tissue in hpatissues]


def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


@follows(mkdir("annotations.dir"))
@split(dbname, ["annotations.dir/*%s" % r for r in outfilesuffixes])
def getDBAnnotations(infile, outfiles):
    '''
    Takes a database (generated using get_gene_annotations.py) and uses this
    to build a series of AnnotationSets.

    One AnnotationSet is generated for each table in the database with
    the $annot suffix.

    AnnotationSets are stored in output files in the annotations.dir
    directory.
    '''
    dbname = PARAMS['db_name']
    PipelineEnrichment.getDBAnnotations(infile, outfiles, dbname, submit=True)


@follows(getDBAnnotations)
@originate(unmappedouts)
def mapUnmappedAnnotations(outfiles):
    '''
    Allows the user to easily add annotations not in the database.
    Requires an "annotations2annotations.py" command specified in the
    pipeline.ini providing details of a flat file containing the necessary
    information to build the AnnotationSet.
    '''
    ua = "annotations.dir/"
    outstem = outfiles[0].replace(ua, "")
    outstem = outstem.replace("_genestoterms.tsv", "")
    substatement = unmapped[outstem]
    outstem2 = outfiles[0].replace("_genestoterms.tsv", "")
    PipelineEnrichment.getFlatFileAnnotations(substatement, outstem2, dbname)


@follows(mkdir("clean_foregrounds.dir"))
@transform("foregrounds.dir/*.tsv", regex("foregrounds.dir/(.*).tsv"),
           r"clean_foregrounds.dir/\1.tsv")
def cleanForegrounds(infile, outfile):
    '''
    Removes duplicates from the foreground and converts IDs to ensemblg.
    '''
    idtype = PARAMS['id_type']
    dbname = PARAMS['db_name']
    E.info(idtype)
    PipelineEnrichment.cleanGeneLists(infile, outfile, idtype, dbname,
                                      submit=True)


@follows(mkdir("clean_backgrounds.dir"))
@transform("backgrounds.dir/*.tsv",
           regex("backgrounds.dir/(.*).tsv"),
           r"clean_backgrounds.dir/\1.tsv")
def cleanUserBackgrounds(infile, outfile):
    '''
    Removes duplicates from user specified backgrounds and converts IDs to
    ensemblg.
    '''
    idtype = PARAMS['id_type']
    dbname = PARAMS['db_name']
    PipelineEnrichment.cleanGeneLists(infile, outfile, idtype, dbname,
                                      submit=True)


@follows(cleanUserBackgrounds)
@active_if(int(PARAMS['hpa_run']) == 1)
@originate(hpatissues)
def buildHPABackground(outfile):
    '''
    Builds a background geneset based on human protein atlas expression values
    specified in pipeline.ini - allows the user to use a tissue specific
    background
    '''
    tissue = outfile.split("/")[1].split("_")[0]
    PipelineEnrichment.HPABackground(tissue,
                                     PARAMS['hpa_minlevel'],
                                     PARAMS['hpa_supportive'],
                                     outfile,
                                     submit=True)


@follows(mapUnmappedAnnotations)
@merge("annotations.dir/*_genestoterms.tsv",
       "clean_backgrounds.dir/allgenes.tsv")
def buildStandardBackground(infiles, outfile):
    '''
    Builds a standard background of all possible genes with any kind of
    annotation.  Genes not mapped for a specific annotation source
    are omitted later, during enrichment testing for that source.
    '''
    statement = "cut -f1 %s | sort | uniq > %s" % (" ".join(infiles), outfile)
    P.run()


@follows(buildHPABackground)
@follows(cleanUserBackgrounds)
@follows(cleanForegrounds)
@follows(mapUnmappedAnnotations)
@follows(buildStandardBackground)
@follows(mkdir("results.dir"))
@product("clean_backgrounds.dir/*",
         formatter(".+/(?P<NAM>.*).tsv"),
         "clean_foregrounds.dir/*",
         formatter(".+/(?P<NAM>.*).tsv"),
         "annotations.dir/*_genestoterms.tsv",
         formatter(".+/(?P<NAM>.*)_genestoterms.tsv"),
         [r'results.dir/{NAM[0][0]}_v_{NAM[1][0]}_{NAM[2][0]}.tsv',
          r'results.dir/{NAM[0][0]}_v_{NAM[1][0]}_{NAM[2][0]}_sig.tsv'])
def foregroundsVsBackgrounds(infiles, outfiles):
    '''
    Takes every possible set of one foreground, one background and one
    AnnotationSet and performs enrichment analysis.  Analysis is
    performed based on the "stats" parameters in the pipeline.ini.
    Results are written to tab delimited files in results.dir, the _sig.tsv
    output file contains signficantly enriched terms only, the other output
    file contains all terms.
    '''
    PipelineEnrichment.foregroundsVsBackgrounds(infiles,
                                                outfiles[0], outfiles[1],
                                                PARAMS['stats_testtype'],
                                                PARAMS['stats_runtype'],
                                                PARAMS['stats_correction'],
                                                PARAMS['stats_thresh'],
                                                dbname,
                                                int(PARAMS[
                                                    'stats_writegenes']),
                                                PARAMS['db_species'],
                                                int(PARAMS['stats_ngenes']),
                                                PARAMS['id_type'],
                                                submit=True)


@follows(mkdir("cytoscape.dir"))
@transform(foregroundsVsBackgrounds,
           regex("results.dir/(.*)(_go|_reactome)(.tsv)"),
           r"cytoscape.dir/\1\2.tsv")
def makeCytoscapeInputs(infiles, outfile):
    infile = infiles[1]
    T = P.getTempFilename(".")
    statement = """
    awk -F "\\t" '{printf("%%%%s\\t%%%%s\\t%%%%s\\t%%%%s\\t+1\\n",\
    $1, $11, $7, $8)}' %(infile)s > %(T)s""" % locals()
    P.run()
    typ = infile.split("_")[-2]
    E.info(typ)
    keep = [line.strip() for line in
            IOTools.openFile(PARAMS['cytoscape_%s' % typ]).readlines()]
    tab = pd.read_csv(T, sep="\t")
    tab = tab[tab['term_id'].isin(keep)]
    tab.columns = ['ID', 'Description', 'pvalue', 'padj', 'Phenotype']
    tab.to_csv(outfile, sep="\t", index=None)
    os.remove(T)


@follows(foregroundsVsBackgrounds)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
