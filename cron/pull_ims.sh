#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

rootdir="${SDPC_ANCILLARY_ROOT}/ims"
subdir="$(date -u +%Y)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
   prev_path=""
else
   prev_file=$(ls -t $target_dir | head -n1)
   prev_path="$target_dir/$prev_file"
fi

lftp ftp://sidads.colorado.edu <<- EOF
   set xfer:log-file
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   cd pub/DATASETS/NOAA/G02156/GIS/1km
   mirror -c --newer-than=now-4days $subdir $target_dir
   quit
EOF

newer=""
if test x"$prev_path" != x ; then
   newer="-newer $prev_path"
fi
new_list=$(find $target_dir -name "ims???????_1km_GIS_v*.tif.gz" $newer)

if test x"$new_list" != x ; then
   export SDPC_ANCILLARY_IMS_DBFILE="$rootdir/ims.sqlite"
   $SDPC_ANCILLARY_ROOT/src/register_ims.py $new_list
fi

