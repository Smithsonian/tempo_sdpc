#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

rootdir="${SDPC_ANCILLARY_ROOT}/goes"
subdir="$(date -u +%Y/%j)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

lftp AO_TEMPO_SERV1@140.90.190.143 <<- EOF
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set ssl:verify-certificate no
   set mirror:require-source true
   set mirror:sort-by name
   set mirror:order *.sha1 *.nc
   mirror -c -O $target_dir -F /PDAFileLinks/g1?_cmi
   quit
EOF

# The 'today' symlink always points to today's GOES data
cd $rootdir
ln -nf -s $subdir today
