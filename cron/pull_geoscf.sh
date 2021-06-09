#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

: "${TESTARGS:=}"

rootdir="${SDPC_ANCILLARY_ROOT}/geos_cf"
subdir="$(date $TESTARGS +%Y/%j)"
daytag="$(date $TESTARGS +%Y%m%d)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

source_url="https://gmao.gsfc.nasa.gov/gmaoftp/kknowlan/TEMPO"

filename_regex="GEOS-CF.v01.rpl.sat_inst_1hr_r720x361_v72.${daytag}_????z.nc4"

lftp <<- EOF
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   mget -c -P 4 -O $target_dir $source_url/$filename_regex
   quit
EOF

# The 'today' symlink always points to today's GEOS-CF data
cd $rootdir
ln -nf -s $subdir today
