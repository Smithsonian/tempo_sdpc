#! /bin/bash

set -e
set -u

if test $# -lt 3 ; then
    echo "Usage: $(basename $0) destdir <source-url> geoscf_version"
    exit 0
fi
rootdir="$1"
source_url="$2"
geoscf_version="$3"

# default to forecast download
want_analysis=0
case "$geoscf_version" in
  *ana* )
     want_analysis=1
     ;;
  *)
     ;;
esac

# To simply tracking what we've downloaded,
# we first download files to an 'incoming' directory.
incoming_dir="$rootdir/incoming"
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir
fi
export incoming_dir

# Any remaining command line args are passed to 'date'.
# For example, to pull data for a particular date, do something like:
#     pull_geoscf.sh <destdir> <source-url> <version> 2020-07-01
today_args=""
today_opt=""
yesterday_opt="-d yesterday"
tomorrow_opt="-d tomorrow"
if test $# -gt 3; then
   shift 3
   today_args="$@"
   today_opt="-d $@"
   yesterday_opt="-d ${today_args}-1day"
   tomorrow_opt="-d ${today_args}+1day"
fi

if test $want_analysis -eq 0 ; then
   # Always get forecasts based on the previous day's replay:
   rpl_day="$(date -u $yesterday_opt +%Y%m%d)"
   source_url="$source_url/$(date -u $yesterday_opt +Y%Y/M%m/D%d)"
else
   # For analysis, we download only today's:
   rpl_day="$(date -u $today_opt +%Y%m%d)"
   source_url="$source_url/$(date -u $today_opt +Y%Y/M%m/D%d)"
fi
export source_url

# Re-order the file variable dimensions
reorder_dims()
{
   path=$1
   # replace '+' with '_' because NCO tools refuses to deal with filenames that contain '+'
   clean_basename="$(basename $path .nc4 | tr + _)"
   infile="$incoming_dir/${clean_basename}.nc4"
   outfile="$incoming_dir/${clean_basename}_reorder.nc4"
   if test x"$path" != x"$infile" ; then
      /bin/mv $path $infile
   fi
   ncpdq -O -a "time,lon,lat,-lev" --no_tmp_fl "$infile" "$outfile" && /bin/rm -f "$infile"
}
export -f reorder_dims

fetch_forecast_for_date()
{
  date_opt="$1"

  fcst_day="$(date -u $date_opt +%Y%m%d)"

  case "$geoscf_version" in
    1)
      fcst_root="GEOS-CF.v01.fcst.sat_inst_1hr_r721x361_v72.${rpl_day}_12z+${fcst_day}"
      fcst_regex="${fcst_root}_??00z.nc4"
      fcst_fmt="H12/${fcst_root}_%sz.nc4"
      ;;

    2)
      fcst_root="GEOS.cf.fcst.sat_inst_1hr_reg_L721x361_v72.${rpl_day}_09z+${fcst_day}"
      fcst_regex="${fcst_root}_??00z.R0.nc4"
      fcst_fmt="${fcst_root}_%sz.R0.nc4"
      ;;

    2R1-ana)
      fcst_root="GEOS.cf.ana.sat_inst_1hr_reg_L721x361_v72.${rpl_day}"
      fcst_regex="${fcst_root}_??00z.R1.nc4"
      fcst_fmt="${fcst_root}_%sz.R1.nc4"
      ;;

    *)
      echo "*** Error: unsupported GEOS-CF version: $geoscf_version"
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
if test $want_analysis -eq 0 ; then
  fetch_forecast_for_date "$tomorrow_opt"
fi
