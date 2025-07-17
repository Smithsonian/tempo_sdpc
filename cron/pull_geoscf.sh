#! /bin/bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_GEOSCF_VERSION:=1}"

set -e
set -u

if test $# -lt 1 ; then
    echo "Usage: $0 <source-url>"
    exit 0
fi
source_url="$1"

rootdir="${SDPC_ANCILLARY_ROOT}/var/geoscf"

# To simply tracking what we've downloaded,
# we first download files to an 'incoming' directory.
incoming_dir="$rootdir/incoming"
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir
fi
export incoming_dir

# Any command line args are passed to 'date'.
# For example, to pull data for a particular date, do something like:
#     pull_geoscf.sh <source-url> 2020-07-01
today_args=""
today_opt=""
yesterday_opt="-d yesterday"
tomorrow_opt="-d tomorrow"
if test $# -gt 1; then
   shift
   today_args="$@"
   today_opt="-d $@"
   yesterday_opt="-d ${today_args}-1day"
   tomorrow_opt="-d ${today_args}+1day"
fi

# Always get forecasts based on the previous day's replay:
rpl_day="$(date -u $yesterday_opt +%Y%m%d)"
source_url="$source_url/$(date -u $yesterday_opt +Y%Y/M%m/D%d)"
export source_url

# Re-order the file variable dimensions
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
export -f reorder_dims

fetch_forecast_for_date()
{
  date_opt="$1"

  fcst_day="$(date -u $date_opt +%Y%m%d)"

  case "$SDPC_GEOSCF_VERSION" in
    1)
      fcst_root="GEOS-CF.v01.fcst.sat_inst_1hr_r721x361_v72.${rpl_day}_12z+${fcst_day}"
      fcst_regex="${fcst_root}_??00z.nc4"
      fcst_fmt="H12/${fcst_root}_%sz.nc4"
      ;;

    2)
      # Sample files provided for testing were initialized at 21z,
      # but production files are expected to be initialized at 09z.
      fcst_root="GEOS.cf.fcst.sat_inst_1hr_reg_L721x361_v72.${rpl_day}_09z+${fcst_day}"
      fcst_regex="${fcst_root}_??00z.V01.nc4"
      fcst_fmt="${fcst_root}_%sz.V01.nc4"
      ;;

    *)
      echo "*** Error: unsupported GEOS-CF version: $SDPC_GEOSCF_VERSION"
      return 1
  esac

  # After reformatting, the files will be moved
  # to their final location:
  subdir="$(date -u $date_opt +%Y/%j)"
  target_dir="$rootdir/$subdir"
  if ! test -d $target_dir ; then
     mkdir -p $target_dir
  fi

  tmpscript=$(mktemp)
  cat <<- EOF > $tmpscript
	set log:file/xfer ""
	set xfer:use-temp-file yes
	set xfer:temp-file-name *.lftp
	set ssl:verify-certificate no
	EOF
  HOURS="$(seq -w 0000 100 2300)"
  for h in $HOURS ; do
      printf "repeat --until-ok -d 5 -c 5 get -c -O $incoming_dir $source_url/$fcst_fmt\n" $h >> $tmpscript
  done
  echo exit >> $tmpscript

  # Download using lftp:
  lftp -f $tmpscript && /bin/rm $tmpscript

  # reorder file variable dimensions
  find $incoming_dir -name $fcst_regex | parallel --will-cite --max-procs 12 reorder_dims {}

  # Move reordered files to their final location, and
  # register them in the asdc upload database
  files=$(find $incoming_dir -name "*_reorder.nc4")
  for f in $files ; do
     final_path="$target_dir/$(basename $f)"
     /bin/mv $f $final_path
     asdc_files.py --dbfile $rootdir/geoscf.sqlite --add $final_path
  done
}

fetch_forecast_for_date "$today_opt"
fetch_forecast_for_date "$tomorrow_opt"
