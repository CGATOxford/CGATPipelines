#!/usr/bin/env bash

# Thanks to:
# https://unix.stackexchange.com/questions/27013/displaying-seconds-as-days-hours-mins-seconds

if [[ $# -ne 1 ]] ; then
   echo
   echo " ERROR: ./secs2hours.sh <seconds>"
   echo
   exit 1
fi

T=$1
D=$((T/60/60/24))
H=$((T/60/60%24))
M=$((T/60%60))
S=$((T%60))
(( $D > 0 )) && printf '%d days ' $D
(( $H > 0 )) && printf '%d hours ' $H
(( $M > 0 )) && printf '%d minutes ' $M
(( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
printf '%d seconds\n' $S
