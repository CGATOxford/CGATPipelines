#!/bin/bash

# Function to print pretty output
pprint() {

   printf "%s\t%s\t%s\n" $1 $2 $3
}

# sanity check; we want input files
if [ $# -eq 0 ] ; then

   echo
   echo " ERROR: no input files given "
   echo " USAGE: merge_testing_output.sh test_mapping.exist test_mapping.lines test_mapping.md5"
   echo
   exit 1;

fi

# check which input files we have
md5file=`echo $@ | sed 's/\s/\n/g' | grep md5`
linesfile=`echo $@ | sed 's/\s/\n/g' | grep lines`
existfile=`echo $@ | sed 's/\s/\n/g' | grep exist`

# check whether they can be read
md5check=0 && linescheck=0 && existcheck=0
[ "$md5file" ] || [ ! -r $md5file ] || [ ! -s $md5file ] && md5check=1
[ "$linesfile" ] || [ ! -r $linesfile ] || [ ! -s $linesfile ] && linescheck=1
[ "$existfile" ] || [ ! -r $existfile ] || [ ! -s $existfile ] && existcheck=1

if [ -z "$md5check" ] && [ -z "$linescheck" ] && [ -z "$existcheck" ] ; then

   echo
   echo " ERROR: one or more input files do not exist, are empty or cannot be read"
   echo " Input files were: " $@
   echo
   exit 1;

fi

# merge all file names
awk '{print $1;}' $md5file $linesfile $existfile | sort -u > temp

pprint "filename" "nlines" "md5"

for file in `cat temp` ; do

   awkfile=$(echo $file | sed 's/\//\\\//g')
   [ -n "$md5file"   ] && md5=`awk '/'$awkfile'/ { print $2; }' $md5file`
   [ -n "$linesfile" ] && nlines=`awk '/'$awkfile'/ { print $2; }' $linesfile`

   [ -z "$nlines" ] && [ -z "$md5" ] && pprint $file "-1" "-1"
   [ -z "$nlines" ] && [ -n "$md5" ] && pprint $file "-1" $md5
   [ -n "$nlines" ] && [ -n "$md5" ] && pprint $file $nlines $md5
   [ -n "$nlines" ] && [ -z "$md5" ] && pprint $file $nlines "-1"

done

# cleanup temp file
rm -f temp

