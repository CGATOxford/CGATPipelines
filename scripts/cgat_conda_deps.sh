#!/usr/bin/env bash
#
# Script to find Python and R deps in this repository
#
# It uses snakefood to find the python dependencies
#

# exit when a command fails
#set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
#set -o nounset

# trace what gets executed
#set -o xtrace


# Global variables

SCRIPT_NAME="$0"
SCRIPT_OPTS="$@"
SCRIPT_FOLDER=$(dirname $0)
REPO_FOLDER=$(dirname ${SCRIPT_FOLDER})
TMP_SFOOD=$(mktemp)
TMP_GREP=$(mktemp)
TMP_MISC=$(mktemp)
TMP_DEPS=$(mktemp)

# dictionary with production CGAT pipelines
declare -A CGAT_PIPELINES
CGAT_PIPELINES[annotations]="pipeline_annotations.py"
CGAT_PIPELINES[enrichment]="pipeline_enrichment.py"
CGAT_PIPELINES[intervals]="pipeline_intervals.py"
CGAT_PIPELINES[mapping]="pipeline_mapping.py"
CGAT_PIPELINES[peakcalling]="pipeline_peakcalling.py"
CGAT_PIPELINES[readqc]="pipeline_readqc.py"
CGAT_PIPELINES[rnaseqde]="pipeline_rnaseqdiffexpression.py"
CGAT_PIPELINES[rnaseqqc]="pipeline_rnaseqqc.py"
CGAT_PIPELINES[scrnaseqqc]="pipeline_scrnaseqqc.py"
CGAT_PIPELINES[windows]="pipeline_windows.py"

# dictionary to translate Python deps
declare -A PY_DEPS
PY_DEPS[Bio]="biopython"
PY_DEPS[MySQLdb]="mysqlclient"
PY_DEPS[SphinxReport]="ignore"
PY_DEPS[alignlib_lite]="alignlib-lite"
PY_DEPS[brewer2mpl]="brewer2mpl"
PY_DEPS[bs4]="beautifulsoup4"
PY_DEPS[bx]="bx-python"
PY_DEPS[configparser]="ignore"
PY_DEPS[drmaa]="python-drmaa"
PY_DEPS[future]="future"
PY_DEPS[ggplot]="ggplot"
PY_DEPS[httplib2]="httplib2"
PY_DEPS[intermine]="intermine"
PY_DEPS[jinja2]="jinja2"
PY_DEPS[lzo]="python-lzo"
PY_DEPS[matplotlib]="matplotlib"
PY_DEPS[metaphlan_utils]="ignore"
PY_DEPS[mygene]="mygene"
PY_DEPS[networkx]="networkx"
PY_DEPS[numpy]="numpy"
PY_DEPS[openpyxl]="openpyxl"
PY_DEPS[pandas]="pandas"
PY_DEPS[pika]="pika"
PY_DEPS[psycopg2]="psycopg2"
PY_DEPS[pyBigWig]="pybigwig"
PY_DEPS[pybedtools]="pybedtools"
PY_DEPS[pysam]="pysam"
PY_DEPS[pyximport]="ignore"
PY_DEPS[rdflib]="rdflib"
PY_DEPS[rpy2]="rpy2"
PY_DEPS[ruffus]="ruffus"
PY_DEPS[scipy]="scipy"
PY_DEPS[seaborn]="seaborn"
PY_DEPS[six]="six"
PY_DEPS[sklearn]="scikit-learn"
PY_DEPS[sqlalchemy]="sqlalchemy"
PY_DEPS[toposort]="toposort"
PY_DEPS[web]="web.py"
PY_DEPS[weblogolib]="python-weblogo"
PY_DEPS[bashlex]="ignore"


# dictionary to translate R deps
declare -A R_DEPS
R_DEPS[BiSeq]="bioconductor-biseq"
R_DEPS[Biobase]="bioconductor-biobase"
R_DEPS[ChIPQC]="bioconductor-chipqc"
R_DEPS[DESeq2]="bioconductor-deseq2"
R_DEPS[DESeq]="bioconductor-deseq"
R_DEPS[DEXSeq]="bioconductor-dexseq"
R_DEPS[GMD]="r-gmd"
R_DEPS[HiddenMarkov]="r-hiddenmarkov"
R_DEPS[HilbertVis]="bioconductor-hilbertvis"
R_DEPS[Hmisc]="r-hmisc"
R_DEPS[IHW]="bioconductor-ihw"
R_DEPS[KEGGdb]="bioconductor-kegg.db"
R_DEPS[MASS]="r-mass"
R_DEPS[MEDIPS]="bioconductor-medips"
R_DEPS[RColorBrewer]="r-rcolorbrewer"
R_DEPS[RSQLite]="r-rsqlite"
R_DEPS[VennDiagram]="r-venndiagram"
R_DEPS[WGCNA]="r-wgcna"
R_DEPS[affy]="bioconductor-affy"
R_DEPS[amap]="r-amap"
R_DEPS[biomaRt]="bioconductor-biomart"
R_DEPS[coloc]="r-coloc"
R_DEPS[cummeRbund]="bioconductor-cummerbund"
R_DEPS[database]="ignore"
R_DEPS[dplyr]="r-dplyr"
R_DEPS[edgeR]="bioconductor-edger"
R_DEPS[flashClust]="r-flashclust"
R_DEPS[gcrma]="bioconductor-gcrma"
R_DEPS[genome_file]="ignore"
R_DEPS[ggplot2]="r-ggplot2"
R_DEPS[gplots]="r-gplots"
R_DEPS[gridExtra]="r-gridextra"
R_DEPS[grid]="r-gridbase"
R_DEPS[gtools]="r-gtools"
R_DEPS[hpar]="bioconductor-hpar"
R_DEPS[limma]="bioconductor-limma"
R_DEPS[maSigPro]="bioconductor-masigpro"
R_DEPS[mapdata]="r-mapdata"
R_DEPS[maps]="r-maps"
R_DEPS[metagenomeSeq]="bioconductor-metagenomeseq"
R_DEPS[optparse]="r-optparse"
R_DEPS[plotrix]="r-plotrix"
R_DEPS[plyr]="r-plyr"
R_DEPS[pvclust]="r-pvclust"
R_DEPS[qqman]="r-qqman"
R_DEPS[reshape2]="r-reshape2"
R_DEPS[reshape]="r-reshape"
R_DEPS[rtracklayer]="bioconductor-rtracklayer"
R_DEPS[samr]="r-samr"
R_DEPS[scales]="r-scales"
R_DEPS[sciplot]="r-sciplot"
R_DEPS[siggenes]="bioconductor-siggenes"
R_DEPS[simpleaffy]="bioconductor-simpleaffy"
R_DEPS[sleuth]="r-sleuth"
R_DEPS[snow]="r-snow"
R_DEPS[spp]="ignore"
R_DEPS[vegan]="r-vegan"
R_DEPS[wasabi]="r-wasabi"
R_DEPS[zinba]="ignore"


# dictionary to translate Misc deps
declare -A MISC_DEPS
MISC_DEPS[AddOrReplaceReadGroups]="picard"
MISC_DEPS[BioProspector]="ignore" # ignore: https://github.com/CGATOxford/CGATPipelines/issues/355
MISC_DEPS[BroadPeak]="ignore" # ignore: https://github.com/CGATOxford/CGATPipelines/issues/355
MISC_DEPS[CalculateHsMetrics]="picard"
MISC_DEPS[CollectAlignmentSummaryMetrics]="picard"
MISC_DEPS[CollectGcBiasMetrics]="picard"
MISC_DEPS[CollectInsertSizeMetrics]="picard"
MISC_DEPS[CollectMultipleMetrics]="picard"
MISC_DEPS[CollectRnaSeqMetrics]="picard"
MISC_DEPS[GenomeAnalysisTK]="gatk"
MISC_DEPS[MarkDuplicates]="picard"
MISC_DEPS[MergeSamFiles]="picard"
MISC_DEPS[PhyloCSF]="ignore" # pipeline_rnaseqlncrna.py is not in production (missing tests)
MISC_DEPS[SICER-rb.sh]="sicer"
MISC_DEPS[SICER.sh]="sicer"
MISC_DEPS[STAR]="star"
MISC_DEPS[ShortStack]="shortstack"
MISC_DEPS[SnpSift.sh]="ignore" # pipeline_exome.py is not in production (missing tests)
MISC_DEPS[ascp]="ignore"
MISC_DEPS[axtToMaf]="ignore" # pipeline_rnaseqlncrna.py is not in production (missing tests)
MISC_DEPS[bamToBed]="bedtools"
MISC_DEPS[bcftools]="ignore" # pipeline_exome.py is not in production (missing tests)
MISC_DEPS[bedGraphToBigWig]="ucsc-bedgraphtobigwig"
MISC_DEPS[bedToBigBed]="ucsc-bedtobigbed"
MISC_DEPS[bedtools]="bedtools"
MISC_DEPS[bgzip]="htslib"
MISC_DEPS[bismark]="bismark"
MISC_DEPS[bismark_methylation_extractor]="ignore" # pipeline_rrbs.py is not in production (missing tests)
MISC_DEPS[bits_test]="ignore"
MISC_DEPS[bowtie-build]="bowtie"
MISC_DEPS[bowtie2-build]="bowtie2"
MISC_DEPS[bowtie2]="bowtie2"
MISC_DEPS[bowtie]="bowtie"
MISC_DEPS[butter]="ignore" # butter is no longer used
MISC_DEPS[bwa]="bwa"
MISC_DEPS[calc]="ignore" # pipeline_exome.py is not in production (missing tests)
MISC_DEPS[cassi]="ignore" # pipeline_gwas.py is not in production (missing tests)
MISC_DEPS[cat]="coreutils"
MISC_DEPS[complementBed]="bedtools"
MISC_DEPS[cuffcompare]="cufflinks"
MISC_DEPS[cuffdiff]="cufflinks"
MISC_DEPS[cufflinks]="cufflinks"
MISC_DEPS[cuffmerge]="cufflinks"
MISC_DEPS[diamond]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[fasta-get-markov]="ignore" # pipeline_motifs.py is not in production (missing tests)
MISC_DEPS[fastq-to-fasta.py]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[fastq_screen]="fastq-screen"
MISC_DEPS[fastqc]="fastqc"
MISC_DEPS[fastx_collapser]="fastx_toolkit"
MISC_DEPS[fastx_reverse_complement]="fastx_toolkit"
MISC_DEPS[featureCounts]="subread"
MISC_DEPS[flash]="flash"
MISC_DEPS[flashpca]="ignore" # pipeline_gwas.py is not in production (missing tests)
MISC_DEPS[gat-run.py]="gat"
MISC_DEPS[gemini]="ignore" # pipeline_exome.py is not in production (missing tests)
MISC_DEPS[go2fmt.pl]="ignore" # ignore: https://github.com/CGATOxford/CGATPipelines/issues/355
MISC_DEPS[grep]="grep"
MISC_DEPS[gsnap]="gmap"
MISC_DEPS[gtfToGenePred]="ucsc-gtftogenepred"
MISC_DEPS[gtf_splicesites]="gmap"
MISC_DEPS[gtf_to_fasta]="tophat"
MISC_DEPS[gzip]="ignore" # gzip not available in conda
MISC_DEPS[hisat2]="hisat2"
MISC_DEPS[idr]="idr"
MISC_DEPS[iit_store]="gmap"
MISC_DEPS[intersectBed]="bedtools"
MISC_DEPS[java]="ignore"
MISC_DEPS[jupyter]="jupyter"
MISC_DEPS[kallisto]="kallisto"
MISC_DEPS[kraken-mpa-report]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[kraken]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[lcamapper.sh]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[locuszoom]="ignore" # pipeline_gwas.py is not in production (missing tests)
MISC_DEPS[macs2]="macs2"
MISC_DEPS[macs]="ignore" # in favour of macs2
MISC_DEPS[map2slim]="ignore" # ignore: https://github.com/CGATOxford/CGATPipelines/issues/355
MISC_DEPS[mast]="meme"
MISC_DEPS[meme-chip]="ignore" # pipeline_motifs.py is not in production (missing tests)
MISC_DEPS[meme]="meme"
MISC_DEPS[mergeBed]="bedtools"
MISC_DEPS[metaphlan.py]="ignore" # pipeline_metagenomecommunities.py is not in production
MISC_DEPS[pandaseq]="pandaseq"
MISC_DEPS[paste]="coreutils"
MISC_DEPS[peakranger]="peakranger"
MISC_DEPS[plink2]="ignore" # pipeline_gwas.py is not in production (missing tests)
MISC_DEPS[prefetch]="ignore"
MISC_DEPS[python]="ignore"
MISC_DEPS[rMATS]="ignore" # pipeline_splicing.py is not in production (missing tests)
MISC_DEPS[rm]="coreutils"
MISC_DEPS[rmats2sashimiplot]="ignore" # pipeline_splicing.py is not in production (missing tests)
MISC_DEPS[sailfish]="sailfish"
MISC_DEPS[salmon]="salmon"
MISC_DEPS[samtools]="samtools"
MISC_DEPS[sickle]="sickle-trim"
MISC_DEPS[slopBed]="bedtools"
MISC_DEPS[snpEff.sh]="ignore" # pipeline_exome.py is not in production (missing tests)
MISC_DEPS[solid2fastq]="ignore" # SOLiD sequencing technology is no longer in use
MISC_DEPS[sortBed]="bedtools"
MISC_DEPS[sort]="coreutils"
MISC_DEPS[sorted_bdg]="ignore"
MISC_DEPS[stampy.py]="ignore" # not available in conda
MISC_DEPS[stringtie]="stringtie"
MISC_DEPS[sum_reads]="ignore" # pipeline_rrbs.py is not in production (missing tests)
MISC_DEPS[tabix]="htslib"
MISC_DEPS[tomtom]="meme"
MISC_DEPS[tophat2]="tophat"
MISC_DEPS[tophat]="tophat"
MISC_DEPS[tr]="coreutils"
MISC_DEPS[transfac2meme]="ignore" # pipeline_motifs.py is not in production (missing tests)
MISC_DEPS[trim_galore]="trim-galore"
MISC_DEPS[trimmomatic]="trimmomatic"
MISC_DEPS[vcf-compare]="vcftools"
MISC_DEPS[vcf-isec]="ignore" # pipeline_exome_cancer.py not in production (missing tests)
MISC_DEPS[vcf-stats]="vcftools"
MISC_DEPS[wc]="coreutils"
MISC_DEPS[wget]="wget"
MISC_DEPS[wigToBigWig]="ucsc-wigtobigwig"
MISC_DEPS[zcat]="ignore" # gzip not available in conda


# function to print the name of production pipelines
# http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
print_pipeline_names() {
   for name in "${!CGAT_PIPELINES[@]}";
   do
      echo " - "$name ;
   done
}


# function to report issues and exit
report_problem() {
   echo
   echo $1
   echo 
   echo " Aborting. "
   exit 1
}


# function to find python imports
# output will go to TMP_SFOOD
find_python_imports() {

sfood ${REPO_FOLDER}/$1 2>&1 \
 | grep 'WARNING     :   ' \
 | grep -v Line \
 | awk '{print $3;}' \
 | grep -v '^_.*' \
 | sed 's/\..*//g' \
 | egrep -v 'Bed|CSV|Database|Experiment|FastaIterator|Fastq|GFF|GTF|Histogram|IOTools|IndexedFasta|Nucmer|Stats' \
 | egrep -v 'Pipeline|PipelineGenomeAssembly|PipelineMapping|PipelineMappingQC|PipelineTracks' \
 | egrep -v 'CGAT|RnaseqDiffExpressionReport|titrationReport|trackers' \
 | sort -u \
 >> ${TMP_SFOOD}

}


# function to find python imports
# output will go to TMP_SFOOD
find_python_imports_recursive() {

sfood -f ${REPO_FOLDER}/CGATPipelines/$1 2>&1 \
 | grep 'WARNING     :   ' \
 | grep -v Line \
 | awk '{print $3;}' \
 | grep -v '^_.*' \
 | sed 's/\..*//g' \
 | egrep -v 'Bed|CSV|Database|Experiment|FastaIterator|Fastq|GFF|GTF|Histogram|IOTools|IndexedFasta|Nucmer|Stats' \
 | egrep -v 'Pipeline|PipelineGenomeAssembly|PipelineMapping|PipelineMappingQC|PipelineTracks' \
 | egrep -v 'CGAT|RnaseqDiffExpressionReport|titrationReport|trackers' \
 | egrep -v '_|EasyDialogs|SOCKS|ce|msvcrt|nt|org|os2|riscos|rourl2path' \
 | sort -u \
 >> ${TMP_SFOOD}

}


# function to find R imports
# output will go to TMP_GREP
find_r_imports() {

grep -i 'library(' -r ${REPO_FOLDER}/$1 \
 | egrep -v 'Binary|zinba' \
 | sed -e 's/\(.*\)library\(.*\)/\2/' \
 | sed 's/[()"&,.%'\'']//g' \
 | sed 's/\\n$//g' \
 | egrep '^[a-zA-Z]{2,}' \
 | sort -u \
 >> ${TMP_GREP}

}


# function to find R imports
# output will go to TMP_GREP
find_r_imports_recursive() {

   TMP_EXT=$(mktemp)
   sfood -q ${REPO_FOLDER}/CGATPipelines/$1 \
    | awk '{print $3" "$4}' \
    | awk '/CGATPipelines/ {print $2}' \
    | sed 's/CGATPipelines\///g' \
    | sed 's/)//g' \
    | sed 's/'\''//g' \
    > ${TMP_EXT}

   echo $1 >> ${TMP_EXT}

   for pkg in `cat ${TMP_EXT}` ;
   do

      find_r_imports "CGATPipelines/$pkg"

   done

   # clean up tmp file
   rm ${TMP_EXT}

}


# function to find misc programs
# output will go to TMP_MISC
find_misc_programs() {

   # will use specific py3 env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate py3-basic

   TMP_EXT=$(mktemp)
   find ${REPO_FOLDER}/CGATPipelines -iname "*pipeline*.py" \
        -not -regex '.*pipeline_docs.*' \
        -not -regex '.*template.*' \
        > ${TMP_EXT}

   for code in `cat ${TMP_EXT}` ;
   do

      python ${REPO_FOLDER}/scripts/cgat_check_deps.py ${code} \
       | egrep -v 'PATH|^$|^cgat$|^No|^R|^Rscript|^cd' \
       >> ${TMP_MISC}

   done

   # return unique names
   cat ${TMP_MISC} | sort -u > ${TMP_EXT}
   cp ${TMP_EXT} ${TMP_MISC}

   # revert to original env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate snakefood

   # clean up tmp file
   rm ${TMP_EXT}

}


# function to find misc programs
# output will go to TMP_MISC
find_misc_recursive() {

   TMP_EXT=$(mktemp)
   sfood -q ${REPO_FOLDER}/CGATPipelines/$1 \
    | awk '{print $3" "$4}' \
    | awk '/CGATPipelines/ {print $2}' \
    | sed 's/)//g' \
    | sed 's/'\''//g' \
    | egrep -v '^CGATPipelines\/Pipeline$' \
    > ${TMP_EXT}

   echo "CGATPipelines/$1" >> ${TMP_EXT}

   # will use specific py3 env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate py3-basic

   for code in `cat ${TMP_EXT}` ;
   do
      python ${REPO_FOLDER}/scripts/cgat_check_deps.py ${REPO_FOLDER}/${code} \
       | egrep -v 'PATH|^$|^cgat$|^No|^R|^Rscript|^cd' \
       >> ${TMP_MISC}

   done

   # return unique names
   cat ${TMP_MISC} | sort -u > ${TMP_EXT}
   cp ${TMP_EXT} ${TMP_MISC}

   # revert to original env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate snakefood

   # clean up tmp file
   rm ${TMP_EXT}

}


# function to display help message
help_message() {
   echo
   echo " Scans this repository to look for conda dependencies."
   echo
   echo " To get the dependencies for all pipelines, run:"
   echo " ./cgat_conda_deps.sh --all"
   echo
   echo " To get the dependencies for a specific pipeline, run:"
   echo " ./cgat_conda_deps.sh --pipeline <name>"
   echo
   echo " where <name> is one of the production pipelines: "
   echo
   print_pipeline_names | sort -u
   exit 1
}


# the script starts here

if [[ $# -eq 0 ]] ; then

   help_message

fi

# variable to store the pipeline to scan
PIPELINE_NAME="default"

# parse input parameters

while [[ $# -gt 0 ]]
do
key="$1"

case $key in

    --help)
    help_message
    ;;

    --all)
    PIPELINE_NAME="all"
    shift
    ;;

    --pipeline)
    PIPELINE_NAME="$2"
    shift 2
    ;;

    *)
    help_message
    ;;

esac
done


# requirement: snakefood
source /ifs/apps/conda-envs/bin/activate snakefood

# initialize temp files
cat /dev/null > ${TMP_SFOOD}
cat /dev/null > ${TMP_GREP}
cat /dev/null > ${TMP_MISC}
cat /dev/null > ${TMP_DEPS}


# find deps depending on given input
if [[ "${PIPELINE_NAME}" == "all" ]] || [[ "${PIPELINE_NAME}" == "default" ]] ; then

   # Python
   find_python_imports "CGATPipelines"
   find_python_imports "scripts"

   # R
   find_r_imports "CGATPipelines"
   find_r_imports "scripts"
   find_r_imports "R"

   # Misc
   find_misc_programs

elif [[ ${CGAT_PIPELINES[${PIPELINE_NAME}]+_} ]] ; then
   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   # This means that PIPELINE_NAME exists in the CGAT_PIPELINES array

   # Python
   find_python_imports_recursive ${CGAT_PIPELINES[${PIPELINE_NAME}]}

   # R
   find_r_imports_recursive ${CGAT_PIPELINES[${PIPELINE_NAME}]}

   # Misc
   find_misc_recursive ${CGAT_PIPELINES[${PIPELINE_NAME}]}

else

   echo
   echo " Unknown pipeline to scan: ${PIPELINE_NAME} "
   echo
   echo " Production pipelines are: "
   print_pipeline_names
   report_problem " "

fi


### create header of env file ###

echo
echo "# output generated by ${SCRIPT_NAME} ${SCRIPT_OPTS}"
echo "# on `date`"
echo
echo "name: cgat-p"
echo
echo "channels:"
echo "- bioconda"
echo "- conda-forge"
echo "- defaults"
echo
echo "dependencies:"


### process python deps ###

for pkg in `cat ${TMP_SFOOD}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${PY_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${PY_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${PY_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "? "$pkg >> ${TMP_DEPS}
   fi

done

# print Python section
echo "# python dependencies"
echo "- python"

if [[ "${PIPELINE_NAME}" == "all" ]] || [[ "${PIPELINE_NAME}" == "default" ]] ; then

   # Add others manually:
   echo "- cgat-report" >> ${TMP_DEPS}
   echo "- nose" >> ${TMP_DEPS}
   echo "- pep8" >> ${TMP_DEPS}
   echo "- setuptools" >> ${TMP_DEPS}

   # Print them all sorted
   sed 's/^- bx-python/# WARNING: bx-python is Py2 only but "pip install bx-python" works with Py3/g' ${TMP_DEPS} \
    | sort -u

else

   # Print them all sorted
   sort -u ${TMP_DEPS}

fi


### process R deps ###

cat /dev/null > ${TMP_DEPS}

for pkg in `cat ${TMP_GREP}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${R_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${R_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${R_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "? "$pkg >> ${TMP_DEPS}
   fi

done

# print R section
R_DEPS_SIZE=$(stat --printf="%s" ${TMP_DEPS})

if [[ ${R_DEPS_SIZE} -gt 0 ]] ; then

   echo "# R dependencies"
   echo "- r-base"
   sort -u ${TMP_DEPS} | grep '^\- r'

   # show bioconductor deps
   sort -u ${TMP_DEPS} | grep '^\- bioconductor'

   # show missing deps
   sort -u ${TMP_DEPS} | grep '^? '

fi


### process misc programs ###

cat /dev/null > ${TMP_DEPS}

for pkg in `cat ${TMP_MISC}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${MISC_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${MISC_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${MISC_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "? "$pkg >> ${TMP_DEPS}
   fi

done

# print misc section

# Add these manually, as they are required but don't use the 'statement' variable to run them
echo "- nomkl" >> ${TMP_DEPS}
echo "- zlib" >> ${TMP_DEPS}

echo "# Misc dependencies"
egrep -v 'CGATparameter|BEGIN|END|local_tmpdir|outfile.gz|\||\+|\(|\)|\[|\]|\/|\=|\{|\}|-p|-l|\\t' ${TMP_DEPS} \
 | sed 's/^- macs2/# WARNING: macs2 is Py2 only. Please install it on a separate conda env/g' \
 | sed 's/^- tophat/# WARNING: tophat is Py2 only. Please install it on a separate conda env/g' \
 | sed 's/^- sicer/# WARNING: sicer is Py2 only. Please install it on a separate conda env/g' \
 | sort -u


# Remove temp files
rm ${TMP_SFOOD}
rm ${TMP_GREP}
rm ${TMP_MISC}
rm ${TMP_DEPS}

