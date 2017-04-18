#!/usr/bin/env bash

# References
# http://kvz.io/blog/2013/11/21/bash-best-practices/
# http://jvns.ca/blog/2017/03/26/bash-quirks/

# exit when a command fails
set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
set -o nounset

# trace what gets executed
#set -o xtrace

# Helper function to debug this script
ddebug() {
   if [[ $# -eq 0 ]] ; then
      echo
      echo " ERROR: no input files given for debugging. "
      echo " The ddebug() function needs one input filename to proceed. "
      exit 1;
   fi

   # split filename and extension, with thanks to:
   # http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
   filename=$(basename $1)
   extension="${filename##*.}"
   filename="${filename%.*}"
   outfile=$filename.debug

   # get file sizes
   md5size=`stat --printf="%s" $filename.md5`
   linessize=`stat --printf="%s" $filename.lines`
   existsize=`stat --printf="%s" $filename.exist`

   # save debugging info to .debug file
   echo -e "\nhostname: "`hostname`"\n" > $outfile
   echo -e "md5file: "$md5file >> $outfile
   echo -e "md5size: "$md5size"\n" >> $outfile
   echo -e "linesfile: "$linesfile >> $outfile
   echo -e "linessize: "$linessize"\n" >> $outfile
   echo -e "existfile: "$existfile >> $outfile
   echo -e "existsize: "$existsize"\n" >> $outfile
}

# Function to print pretty output
pprint() {

   printf "%s\t%s\t%s\n" $1 $2 $3
}

# sanity check; we want input files
if [[ $# -eq 0 ]] ; then
   echo
   echo " ERROR: no input files given "
   echo " USAGE: merge_testing_output.sh test_mapping.exist test_mapping.lines test_mapping.md5"
   echo
   exit 1;
fi

# check which input files we have
md5file=`echo $@ | sed 's/\s/\n/g' | grep md5` || md5file=""
linesfile=`echo $@ | sed 's/\s/\n/g' | grep lines` || linesfile=""
existfile=`echo $@ | sed 's/\s/\n/g' | grep exist` || existfile=""

if [[ ! -r "$md5file" ]] && [[ ! -r "$linesfile" ]] && [[ ! -r "$existfile" ]] ; then

   echo
   echo " ERROR: one or more input files do not exist, are empty or cannot be read"
   echo " Input files were: " $@
   echo
   exit 1;

fi

pprint "filename" "nlines" "md5"

# merge all input files and iterate over them
for file in `awk '{print $1;}' $md5file $linesfile $existfile | sort -u` ; do

   awkfile=$(echo $file | sed 's/\//\\\//g')
   md5="" && nlines=""
   [[ -r "$md5file"   ]] && md5=`awk '/'$awkfile'/ { print $2; }' $md5file`
   [[ -r "$linesfile" ]] && nlines=`awk '/'$awkfile'/ { print $2; }' $linesfile`

   [[ -z "$nlines" ]] && [[ -z "$md5" ]] && pprint $file "-1" "-1"
   [[ -z "$nlines" ]] && [[ -n "$md5" ]] && pprint $file "-1" $md5
   [[ -n "$nlines" ]] && [[ -n "$md5" ]] && pprint $file $nlines $md5
   [[ -n "$nlines" ]] && [[ -z "$md5" ]] && pprint $file $nlines "-1"

done

# print debugging info
#ddebug $1

exit 0

