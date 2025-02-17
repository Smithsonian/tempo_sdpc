#! /bin/sh

# When possible, perform  daily INR QA check
# and archive the resulting output.

# Silent exit when the OTS INR QA check script is not available
if ! command -v tempo_inr_quality.sh &> /dev/null ; then
   exit 0
fi
# Silent exit when database doesn't exist
if ! test -f $SDPC_ARCHIVE_DBFILE ; then
   exit 0
fi
# Silent exit when RAD_L1 table doesn't exist
sql_table_query="select name from sqlite_master where type='table' and name='RAD_L1'"
table_name=$(sqlite3 -readonly -cmd ".timeout 20000" $SDPC_ARCHIVE_DBFILE "$sql_table_query")
if test -z "$table_name" ; then
   exit 0
fi

set -u

PGMNAME=$(basename $0)

test_date_ymd=""
if test $# -eq 1 ; then
   test_date_ymd="$1"
fi

root_work_dir="$SDPC_PIPE_DIR/inr/quality"
if ! test -d $root_work_dir ; then
   mkdir -p $root_work_dir || exit 1
fi
# The date most recently processed by INR is stored in this file
date_file="$root_work_dir/date"

target_time="4am"
# The script should be run at least once per day, just after $target_time.
# If it runs more or less often, that will also work, but it's
# designed to generate output around $target_time every day, so it
# does its work only when > 24 hours have elapsed since the last check,
# and it always resets the saved timestamp to aim for that schedule.

tar_prefix="tempo_inrq_d"
# Prefix used to construct tar file name to be archived

tstamp_file="$root_work_dir/timestamp"
# The timestamp file always contains the time_t value for
# $target_time on the day the INR QA script was last run.

# Exit if the INR QA check has been run within the past 24h.
if test -f $tstamp_file ; then
   timet_prev=$(cat $tstamp_file)
   timet_now=$(date +%s)
   if test $(($timet_now - $timet_prev)) -lt 86400 ; then
      exit 0
   fi
fi

# The date to be checked is either stored in a file,
# or provided on the command line. Otherwise, do nothing.
if test x"$test_date_ymd" != x ; then
   date_ymd="$test_date_ymd"
elif test -f "$date_file" ; then
   date_ymd="$(cat $date_file)"
else
   exit 0
fi

have_unchecked_data()
{
   _ymd=$1
   # Unix time_t early on the specified date
   _t1=$(date --date "${_ymd}T04:00:00-06:00" +%s)
   # Unix time_t at the TEMPO epoch (1980-01-06T00:00:00Z)
   _tempo_epoch_timet=315964800
   # TEMPO timestamp early on the specified date
   _tempo_t1="$(($_t1-$_tempo_epoch_timet))"
   _tempo_t2="$(($_tempo_t1+16*3600))"
   # Do we have radiance data in the time range of interest?
   sql_query="select count(*) from RAD_L1 where istart between $_tempo_t1 and $_tempo_t2"
   num_files=$(sqlite3 -readonly -cmd ".timeout 20000" $SDPC_ARCHIVE_DBFILE "$sql_query")
   if test $num_files -eq 0 ; then
      exit 0
   fi
   # Satellite-local day number for the specified date
   _sat_day="$((${_tempo_t1}/86400))"
   # Path to tar file containing output for this date
   _archived_tar_path="$SDPC_ARCHIVE_DIR/L1/RAD/D${_sat_day}/inr/${tar_prefix}${_sat_day}.tar"
   # Silent exit if the tar file is already archived
   if test -f $_archived_tar_path ; then
      exit 0
   fi
}
have_unchecked_data $date_ymd

echo "$PGMNAME: performing INR quality check for date: $date_ymd"

# Update the timestamp file.
date --date "$target_time today" +%s > $tstamp_file

# Perform the INR QA check
log_path="$root_work_dir/log.$(date +%s).txt"
inr_qa_check.sh $date_ymd $root_work_dir > $log_path 2>&1

# Examine log and report any failures
num_fail=$(grep -iwc noncompliant $log_path)
if test $num_fail -gt 0 ; then
   num_ok=$(grep -iwc compliant $log_path)
   num_ind=$(grep -iwc indeterminate $log_path)
   echo "$PGMNAME: WARNING: INR quality check: $date_ymd scans: fail:$num_fail ok:$num_ok indeterminate:$num_ind"
fi

tar_path=$(find $root_work_dir -maxdepth 1 -name "${tar_prefix}?????.tar" -mtime -5)
if ! test -f $tar_path ; then
   exit 0
fi
sat_day=$(basename $tar_path .tar | sed -e s,${tar_prefix},,)

# move the log and tar files into the archive
dest_dir="$SDPC_ARCHIVE_DIR/L1/RAD/D${sat_day}/inr/"
if ! test -d $dest_dir ; then
   mkdir -p $dest_dir
fi
/bin/mv $log_path $dest_dir/${tar_prefix}${sat_day}.log
/bin/mv $tar_path $dest_dir
