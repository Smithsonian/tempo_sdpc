#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $0 <source-url> <source-dir>"
    exit 0
fi
source_url="$1"
source_dir="$2"

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

lftp $source_url <<- EOF
   set xfer:log-file
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   cd $source_dir
   mirror -c --newer-than=now-4days $subdir $target_dir
   quit
EOF

newer=""
if test x"$prev_path" != x ; then
   newer="-newer $prev_path"
fi
new_list=$(find $target_dir -name "ims???????_1km_GIS_v*.tif.gz" $newer)

if test x"$new_list" != x ; then
   $SDPC_ANCILLARY_ROOT/bin/register_ims.py --dbfile $rootdir/ims.sqlite $new_list
fi

