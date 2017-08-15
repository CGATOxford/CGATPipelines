'''
PipelineEnrichmentGSEA.py
=============================================

:Tags: Python

Usage
-----
This pipeline is a wrapper of script runGSEA.py (enrichment analysis
                                                by using GSEA. Further
                                                description is provided
                                                below.)
To run this pipeline, one needs to specify required parameteres in
pipeline.ini file (configuration file).
This pipeline entails steps:
-----------
First step: Preprocessing of gene list(expression data set)
----------- Note: 1. Input gene list should be tab delimited file.
                  	a. First line of dataset will be considered as
                  	   header. Suffix of file name should be ".gene.tsv"
                  	b. Gene ids within gene list and gene set should be the same
                   2. Annotations from a Database:(to convert genelists)
                         a. AnnotationSets are predominantly generated from a database using an
                            AnnotationParser method.
                         b. The Database is generated using the pipeline pipeline_geneinfo.py.
                            This database is required to run pipeline_enrichment.
            Input gene list is translated into required id type.
            (Available options are specified in .ini file), sorts
            the gene list on the basis of provided ranking metric.
            It also removes all duplicate ids and generates report.
            A summary of preprocessing steps of the gene list is provided and lists
	    of duplicate gene ids that were discarded is also listed.
            A new gene list file (after preprocessing is created in a folder
            that has the same name as gene list file name. This new file is used
            for further analysis.
------------
Second step: Call runGSEA.py script file for enrichemnt analysis
-----------
This script will perform the enrichment analysis, by using gene set enrichment analysis
(GSEA) and leading edge analysis.
            "Leading edge are defined as genes that are common to multiple
             significantly enriched gene sets  and  coordinately  enriched
             in a phenotypic comparison of interest.They represent a rich
             source of biologically important  genes."
-----
It takes two input files:
      1. Ranked list of genes (Preprocessed Expression data set file,
                               created by first step of pipeline).
      2. Gene set
          - A gene sets file defines one or more gene sets. For each gene
            set,the file contains the gene set name and the list of genes in
            that gene set. A gene sets file is a tab-delimited text file in
            gmx or gmt format. For descriptions and examples of each file
            format please refer to:
            http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

          - The Molecular Signatures Database (MSigDB)
            (http://software.broadinstitute.org/gsea/msigdb/index.jsp)
            is a collection of annotated gene sets, which can be used for gene
            set enrichment analysis.OR you can create your own gene set in gmt
            or gmx format.
      3. Rest of the parameters can be specified in to pipeline.ini configuration
         file. Every parameter is set to deafult value.


This script will summarize the analysis in the following format:
       1. GSEA Statistics
       2. GSEA Report
       3. Leading edge analysis report

example
-------

The way the test is ran:
cgat runGSEA -f "Expression_data.tsv" -g "Gene_set.gmt" -n 10 -d 1 -l 4

Default run conditions:
cgat runGSEA -f "Expression_data.tsv" -g "Gene_set.gmt"

--------------
GSEA Statistics
---------------
It includes following statistics for GSEA(for each phenotype):
          - Enrichment Score (ES)
          - Normalized Enrichment Score (NES)
          - False Discovery Rate (FDR)
          - Nominal P Value

--------------
GSEA reports
--------------
	- Global Statistics and Plots include:
	    a) Enrichment plot,
            b) Three separate bar plots that provide a quick overview of top 20 (this number is user defined)
	       enriched upregulated, downregulated and overall enriched genesets on the basis of their FDR values.
            c) Global distribution of normalized enrichment score
            d) Global distribution of normalized enrichment score with corresponding FDR q values and p values.
        - Reports:
		1 - Enrichment in Phenotype (of up and downregulated genes)
           	  This report provides summary of enrichment analysis of each phenotype.
            	  It includes details of which genesets are up and downregulated and a summary
	    	  of significant enriched gensets on the basis of FDR and p values.)
         	2 - Gene Set Details
                  This report provides summary of preprocessing steps of the genesets provided and
	   	  lists genes sets that were used in the anlysis and which one were discarded due to set thresholds
         	3 - Detailed Enrichment Results
         	  This report provides detail statistics of each geneset(for each phenotype). Three reports are
         	  generated. report for uoregulated genesets, downregulated genesets, and enriched genesets organised
         	  on the basis of their FDR values.

	    By default, enrichment plot for top 20 gene sets will be reported.

----------------------------
Leading edge analysis report
----------------------------
It will report graphs that help you visualize the overlap between the selected leading edge subsets. It also
summarises the analysis in the form of reports. By default top 10 enriched genesets will be used for leading edge analysis.
	- Leading edge plots include:
             a) Heat Map(unclustered)
          	This provides an overview of overlap between leading edge subsets
             b) Heat Map(clustered)
          	This heat map will be generated after hierarchical clustering of leading edge subset. It will
          	show you clustered genes among subsets
       	     c) Set-to-Set Heat Map
         	This plot help you to visualize intensity of overlap between subsets (i.e. the extent of overlap between two genesets)
	     d) Dendogram to illustrate the arrangement of the clusters produced by hierarchical clustering.

        - Reports:
           1- Leading_edge_summary_report: summary of genesets and corresponding enrichment statistics that were used for the leading edge analysis.
           2- Leading edge matrix (gmx) file provides detailed information on leading edge analysis genesets
              (i.e. participating genes in each gene set).
           3- Leading edge (gct,cluster format) files for unclustered and clustered gene set. It is a boolean matrix.
              that can be used as an input into other resources for additional analysis as this is ideal format for cluster representation
              (in GSEA)

For details on the algorithm please refer to
Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550)
                    and
Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273).

=============================================
'''
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
import numpy as np
import itertools as ITL
import csv
import os
import string
from toposort import toposort_flatten


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
    print(tables)
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


def list_Duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [idx for idx, item in enumerate(
        seq) if item in seen or seen_add(item)]


def read_Expression_data(filename):
    '''
     Read gene list file.
    '''
    f = open(filename, "r")
    express_id = []
    lines = list(csv.reader(f, delimiter="\t"))
    lines.pop(0)
    e_id = [item[0] for item in lines]
    value = [item[1] for item in lines]
    value_arr = np.array(value, dtype=np.float)
    express_value = np.zeros((len(value),), dtype=np.float)
    ind_sort = np.argsort(-value_arr)
    c = 0
    for ii in ind_sort:
        express_value[c] = value_arr[ii]
        express_id.append(e_id[ii])
        c = c + 1
    f.close()
    return express_id, express_value


def readDBTable(dbname, tablename):
    '''
    Reads the specified table from the specified database.
    Returns a list of tuples representing each row
    '''
    dbh = sqlite3.connect(dbname)
    c = dbh.cursor()
    statement = "SELECT * FROM %s" % tablename
    E.warn(statement)
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


def translateGenelist(dbname, genelist, idtype, id_conversion):
    '''
    Translates a list of gene names from idtype to id_conversion based
    on the database table. This table needs to exist in the database.
    '''
    if(id_conversion == "ensemblg"):
        trans = pd.DataFrame(
            readDBTable(
                dbname, "%s2%s$geneid" %
                (id_conversion, idtype)))
        trans.columns = getDBColumnNames(
            dbname, "%s2%s$geneid" %
            (id_conversion, idtype))
    else:
        trans = pd.DataFrame(
            readDBTable(
                dbname, "%s2%s$geneid" %
                (idtype, id_conversion)))
        trans.columns = getDBColumnNames(
            dbname, "%s2%s$geneid" %
            (idtype, id_conversion))
    #print(trans.loc[trans[1] == '29924'])
    mergeon = set(trans.columns)
    mergeon.remove(id_conversion)
    mergeon = list(mergeon)[0]
    database_list = trans[mergeon].tolist()
    index_in_db = [database_list.index(x)
                   for x in genelist if x in database_list]
    newgenelist = trans[id_conversion][index_in_db]
    return list(newgenelist.values)


def create_File_after_preprocessing(file_ids, values, index_to_kept, outfile):
    f = open(outfile, "w")
    f.write("Gene Id\tValues\n")
    for i in index_to_kept:
        f.write(file_ids[i] + "\t" + str(values[i]) + "\n")
    f.close()
    return


def generate_Report_of_preprocessing(file_ids, index_to_remove, report_file):
    f = open(report_file, "w")
    if(len(index_to_remove) > 0):
        f.write("There were duplicate row identifiers in the ranked list." +
                "One id was choosen. Details are below:\n# of row ids" +
                " in original dataset:" + str(len(file_ids)) + "\n# of row UNIQUE ids"
                " in original dataset:" + str((len(file_ids) - len(index_to_remove))) +
                "\n The duplicates were:\n")
        for i in index_to_remove:
            f.write(file_ids[i] + "\n")
        f.close()
    else:
        f.write("There were no duplicate row identifiers in the ranked list." +
                "Details are below:\n# of row ids" +
                " in original dataset:" + str(len(file_ids)) + "\n# of row UNIQUE ids"
                " in original dataset:" + str((len(file_ids) - len(index_to_remove))) + "\n")
        f.close()
    return


def preprocess_ExpressionData(
        filename,
        dbname,
        idtype,
        id_conversion,
        outfile):
    '''
    Preprocess expression dataset: Remove duplicates,Translate a list of
    gene names from idtype to id_conversion and Rank them on the basis
    of ranking metrices
    '''
    filename_base = os.path.basename(filename)
    part1, suff1, suff2 = filename_base.split('.')
    report_file = "_".join([part1, "preprocessing_report.txt"])
    file_ids, values = read_Expression_data(filename)
    # print(len(file_ids))
    if(idtype != id_conversion):
        E.warn("It should not come here")
        E.warn(idtype)
        E.warn(id_conversion)
        file_ids = translateGenelist(dbname, file_ids, idtype, id_conversion)
    '''Remove duplicats.'''
    # print(len(file_ids))
    index_to_kept = list(range(len(file_ids)))
    index_to_remove = list_Duplicates(file_ids)
    # print(len(index_to_remove))
    if(len(index_to_remove) > 0):
        index_to_kept = [i for j, i in enumerate(
            index_to_kept) if j not in index_to_remove]
    generate_Report_of_preprocessing(file_ids, index_to_remove, report_file)
    create_File_after_preprocessing(file_ids, values, index_to_kept, outfile)
    return

# preprocess_ExpressionData("Expression_data_test.gene.tsv","/ifs/projects/reshma/DATABASE/human_db_110817","entrez","ensemblg","Expression_data_test")
