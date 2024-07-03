#! /bin/bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $0 <source-url> <source-dir>"
    exit 0
fi
source_url="$1"
source_dir="$2"

rootdir="${SDPC_ANCILLARY_ROOT}/var/ims"
subdir="$(date -u +%Y)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

lftp $source_url <<- EOF
   set log:file/xfer ""
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   set ssl:verify-certificate no
   cd $source_dir
   mirror -c --newer-than=now-3days $subdir $target_dir
   quit
EOF

new_list=$(find $target_dir -name "ims???????_1km_v*.nc.gz")

if test x"$new_list" != x ; then
   kept_list=""
   for gzfile in $new_list ; do
       file=$(basename $gzfile .gz)
       if test -f $target_dir/$file ; then
          /bin/rm $gzfile
       else
          gunzip $gzfile
          kept_list="$target_dir/$file $kept_list"
       fi
   done
   if test x"$kept_list" != x ; then
      asdc_files.py --dbfile $rootdir/ims.sqlite --add $kept_list
   fi
fi

