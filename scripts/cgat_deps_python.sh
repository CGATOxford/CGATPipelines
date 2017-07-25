#!/usr/bin/env bash
#
# Script to find Python deps in this repository
#
# It uses snakefood to find the python dependencies
#

SCRIPT_FOLDER=$(dirname $0)
REPO_FOLDER=$(dirname ${SCRIPT_FOLDER})

# requirement: snakefood
source /ifs/apps/conda-envs/bin/activate snakefood

sfood ${REPO_FOLDER}/{CGATPipelines,scripts} 2>&1 \
 | grep 'WARNING     :   ' \
 | grep -v Line \
 | awk '{print $3;}' \
 | grep -v '^_.*' \
 | sed 's/\..*//g' \
 | egrep -v 'Bed|CSV|Database|Experiment|FastaIterator|Fastq|GFF|GTF|Histogram|IOTools|IndexedFasta|Nucmer|Stats' \
 | egrep -v 'Pipeline|PipelineGenomeAssembly|PipelineMapping|PipelineMappingQC|PipelineTracks' \
 | egrep -v 'CGAT|RnaseqDiffExpressionReport|titrationReport|trackers' \
 | sort -u

