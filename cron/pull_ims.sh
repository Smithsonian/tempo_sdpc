#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

rootdir="${SDPC_ANCILLARY_ROOT}/ims"
subdir="$(date +%Y)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

lftp ftp://sidads.colorado.edu <<- EOF
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   cd pub/DATASETS/NOAA/G02156/GIS/1km
   mirror -c --newer-than=now-4days $subdir $target_dir
   quit
EOF

