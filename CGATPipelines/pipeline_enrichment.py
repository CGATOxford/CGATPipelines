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

"""
from ruffus import *
from ruffus.combinatorics import *
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import sys
import PipelineGSEnrichment_t2 as PipelineEnrichment


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
unmappedouts = [["untranslated_annotations.dir/%s%s" % (u, s)
                 for s in outfilesuffixes]
                for u in unmapped]


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
    '''
    dbname = PARAMS['db_name']
    PipelineEnrichment.getDBAnnotations(infile, outfiles, dbname, submit=True)


@follows(mkdir("untranslated_annotations.dir"))
@originate(unmappedouts)
def mapUnmappedAnnotations(outfiles):
    '''
    '''
    ua = "untranslated_annotations.dir/"
    outstem = outfiles[0].replace(ua, "")
    outstem = outstem.replace("_genestoterms.tsv", "")
    substatement = unmapped[outstem]
    statement = """python %s/annotations2annotations.py
    -m standardise
    %s --dir %s""" % (PARAMS['scriptsdir'], substatement, ua)
    P.run()


@follows(getDBAnnotations)
@follows(mapUnmappedAnnotations)
@collate("untranslated_annotations.dir/*.tsv",
         regex("untranslated_annotations.dir/(.*)_(.*)_(.*).tsv"),
         [r"annotations.dir/\1%s" % x
          for x in outfilesuffixes])
def translateAnnotations(infiles, outfiles):
    '''
    '''
    dbname = PARAMS['db_name']
    PipelineEnrichment.translateAnnotations(infiles, outfiles, dbname,
                                            submit=True)


@follows(mkdir("clean_foregrounds.dir"))
@transform("foregrounds.dir/*.tsv", regex("foregrounds.dir/(.*).tsv"),
           r"clean_foregrounds.dir/\1.tsv")
def cleanForegrounds(infile, outfile):
    '''
    '''
    idtype = PARAMS['foreground_idtype']
    dbname = PARAMS['db_name']
    PipelineEnrichment.cleanGeneLists(infile, outfile, idtype, dbname,
                                      submit=True)


@follows(mkdir("clean_backgrounds.dir"))
@transform("backgrounds.dir/*.tsv",
           regex("backgrounds.dir/(.*).tsv"),
           r"clean_backgrounds.dir/\1.tsv")
def cleanUserBackgrounds(infile, outfile):
    '''
    '''
    idtype = PARAMS['background_idtype']
    dbname = PARAMS['db_name']
    PipelineEnrichment.cleanGeneLists(infile, outfile, idtype, dbname,
                                      submit=True)


@follows(cleanUserBackgrounds)
@active_if(int(PARAMS['hpa_run']) == 1)
@originate("clean_backgrounds.dir/hpa_background.tsv")
def buildHPABackground(outfile):
    '''
    Builds a background geneset based on human protein atlas expression values
    specified in pipeline.ini.
    '''
    PipelineEnrichment.HPABackground(PARAMS['hpa_tissue'],
                                     PARAMS['hpa_minlevel'],
                                     PARAMS['hpa_supportive'],
                                     outfile,
                                     submit=True)


@follows(translateAnnotations)
@merge("annotations.dir/*_genestoterms.tsv",
       "clean_backgrounds.dir/allgenes.tsv")
def buildStandardBackground(infiles, outfile):
    '''
    '''
    statement = "cut -f1 %s | sort | uniq > %s" % (" ".join(infiles), outfile)
    P.run()


@follows(buildHPABackground)
@follows(cleanUserBackgrounds)
@follows(cleanForegrounds)
@follows(translateAnnotations)
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
    '''
    PipelineEnrichment.foregroundsVsBackgrounds(infiles,
                                                outfiles[0], outfiles[1],
                                                PARAMS['stats_testtype'],
                                                PARAMS['stats_runtype'],
                                                PARAMS['stats_correction'],
                                                PARAMS['stats_thresh'],
                                                submit=True)


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
