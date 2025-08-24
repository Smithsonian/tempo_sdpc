#! /usr/bin/env bash

set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $(basename $0) destdir <source-url>"
    exit 0
fi
rootdir="$1"
source_url="$2"

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
   register_iers.py --dbfile $rootdir/iers.sqlite --rename $latest
fi

