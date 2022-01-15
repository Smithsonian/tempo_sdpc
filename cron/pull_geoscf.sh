#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

if test $# -lt 1 ; then
    echo "Usage: $0 <source-url>"
    exit 0
fi
source_url="$1"

# Any command line args are passed to 'date'.
# For example, to pull data for a particular date, do something like:
#     pull_geoscf.sh <source-url> -d 2020-07-01
dateargs=""
if test $# -gt 1; then
   shift
   dateargs="$@"
fi

rootdir="${SDPC_ANCILLARY_ROOT}/var/geoscf"
subdir="$(date -u $dateargs +%Y/%j)"

# forecast files
daytag="$(date -u $dateargs +%Y%m%d)"
source_url="$source_url/$(date -u $dateargs +Y%Y/M%m/D%d/H12)"

# To simply tracking what we've downloaded,
# we first download files to an 'incoming' directory.
incoming_dir="$rootdir/incoming"
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir
fi

# After reformatting, the files will be moved
# to their final location:
target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

# Download using lftp:
filename_regex="GEOS-CF.v01.rpl.sat_inst_1hr_r720x361_v72.${daytag}_????z.nc4"

lftp <<- EOF
   set xfer:log-file
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   mget -c -P 4 -O $incoming_dir $source_url/$filename_regex
   quit
EOF

# Re-order the variable dimensions

reorder_dims()
{
   path=$1
   bn=$(basename $path .nc4)
   new_file="${bn}_reorder.nc4"
   ncpdq -O -a "time,lon,lat,-lev" --no_tmp_fl $path $incoming_dir/$new_file && /bin/rm -f "$path"
}
export incoming_dir
export -f reorder_dims

find $incoming_dir -name $filename_regex | parallel --will-cite --max-procs 12 reorder_dims {}

# Move files to their final location, and
# register them in the asdc upload database

files=$(find $incoming_dir -name "*.nc4")
for f in $files ; do
   bn=$(basename $f)
   final_path="$target_dir/$bn"
   /bin/mv $f $final_path
   asdc_files.py --dbfile $rootdir/geoscf.sqlite --add $final_path
done
