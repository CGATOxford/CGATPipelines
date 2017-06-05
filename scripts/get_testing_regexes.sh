#!/bin/bash
#
# grep all regexes used in the pipeline.ini file of pipeline_testing.py
#
# useful to debug pipeline_testing.py
#

if [ $# -ne 1 ] ; then

   echo
   echo " Usage: get_testing_regexes.sh /pipeline/testing/pipeline.ini"
   echo
   exit 1 ;

fi

egrep '^suffixes|^regex_linecount|^regex_exist' $1 \
 | sed 's/suffixes\=\|regex_linecount\=\|regex_exist\=//g' \
 | sed 's/,/\n/g' | sed 's/''|''/\n/g' \
 | sort -u

