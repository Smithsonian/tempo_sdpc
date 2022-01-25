#! /bin/bash

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
#     pull_geoscf.sh <source-url> 2020-07-01
today_args=""
today_opt=""
yesterday_opt="-d yesterday"
if test $# -gt 1; then
   shift
   today_args="$@"
   today_opt="-d $@"
   yesterday_opt="-d ${today_args}-1day"
fi

rootdir="${SDPC_ANCILLARY_ROOT}/var/geoscf"
subdir="$(date -u $today_opt +%Y/%j)"

# forecast files
source_url="$source_url/$(date -u $yesterday_opt +Y%Y/M%m/D%d/H12)"
rpl_day="$(date -u $yesterday_opt +%Y%m%d)"
fcst_day="$(date -u $today_opt +%Y%m%d)"
fcst_regex="GEOS-CF.v01.fcst.sat_inst_1hr_r721x361_v72.${rpl_day}_12z+${fcst_day}_????z.nc4"

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

lftp <<- EOF
   set log:file/xfer ""
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set mirror:require-source true
   mget -c -P 4 -O $incoming_dir $source_url/$fcst_regex
   quit
EOF

# Re-order the variable dimensions

reorder_dims()
{
   path=$1
   # replace '+' with '_' because NCO tools refuses to deal with filenames that contain '+'
   clean_basename="$(basename $path .nc4 | tr + _)"
   infile="$incoming_dir/${clean_basename}.nc4"
   outfile="$incoming_dir/${clean_basename}_reorder.nc4"
   /bin/mv $path $infile
   ncpdq -O -a "time,lon,lat,-lev" --no_tmp_fl "$infile" "$outfile" && /bin/rm -f "$infile"
}
export incoming_dir
export -f reorder_dims

find $incoming_dir -name $fcst_regex | parallel --will-cite --max-procs 12 reorder_dims {}

# Move files to their final location, and
# register them in the asdc upload database

files=$(find $incoming_dir -name "*.nc4")
for f in $files ; do
   final_path="$target_dir/$(basename $f)"
   /bin/mv $f $final_path
   asdc_files.py --dbfile $rootdir/geoscf.sqlite --add $final_path
done
