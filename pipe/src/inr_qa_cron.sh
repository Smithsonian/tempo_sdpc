#! /bin/sh

# When possible, perform  daily INR QA check
# and archive the resulting output.

PGMNAME=$(basename $0)

# The date to be checked is provided on the command line
date_ymd=""
if test $# -eq 1 ; then
   date_ymd="$1"
else
   echo "Usage: $PGMNAME YYYY-MM-DD"
   exit 1
fi

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

root_work_dir="$SDPC_PIPE_DIR/inr/quality"
if ! test -d $root_work_dir ; then
   mkdir -p $root_work_dir || exit 1
fi

tar_prefix="tempo_inrq_d"
# Prefix used to construct tar file name to be archived

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
   # Exit if the tar file is already archived
   if test -f $_archived_tar_path ; then
      echo "file exists: $_archived_tar_path"
      exit 0
   fi
}
have_unchecked_data $date_ymd

# Perform the INR QA check
log_path="$root_work_dir/log.$(date +%s).txt"
inr_qa_check.sh $date_ymd $root_work_dir > $log_path 2>&1
exit_status="$?"

# Summarize log, warn if necessary
num_ok=0
num_ind=0
num_non=0
if test -f "$log_path" ; then
   num_ok=$(grep -iwc compliant $log_path)
   num_ind=$(grep -iwc indeterminant $log_path)
   num_non=$(grep -iwc noncompliant $log_path)
fi
if test "$exit_status" -ne 0 ; then
   status="ERROR"
elif test $num_non -gt 0 || test $num_ok -lt $num_ind ; then
   status="WARNING"
else
   status="SUCCESS"
fi
echo "$PGMNAME: $status: INR quality check: $date_ymd scans: noncompliant:$num_non compliant:$num_ok indeterminant:$num_ind"

tar_path=$(find $root_work_dir -maxdepth 1 -name "${tar_prefix}?????.tar" -mtime -5)
if ! test -f $tar_path ; then
   exit 0
fi
sat_day=$(basename $tar_path .tar | sed -e s,${tar_prefix},,)

# Move the log and tar files into the archive
dest_dir="$SDPC_ARCHIVE_DIR/L1/RAD/D${sat_day}/inr/"
if ! test -d $dest_dir ; then
   mkdir -p $dest_dir
fi
/bin/mv $log_path $dest_dir/${tar_prefix}${sat_day}.log
/bin/mv $tar_path $dest_dir
