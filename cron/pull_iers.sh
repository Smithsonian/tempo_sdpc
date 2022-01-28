#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

if test $# -ne 1 ; then
    echo "Usage: $0 <source-url>"
    exit 0
fi
source_url="$1"

rootdir="${SDPC_ANCILLARY_ROOT}/var/iers"

target_dir="$rootdir/files"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

latest="$target_dir/bulletina.txt"
lftp <<- EOF
   set log:file/xfer ""
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   get $source_url -o $latest
   quit
EOF

if test x"$latest" != x ; then
   $SDPC_ANCILLARY_ROOT/bin/register_iers.py --dbfile $rootdir/iers.sqlite --rename $latest
fi

