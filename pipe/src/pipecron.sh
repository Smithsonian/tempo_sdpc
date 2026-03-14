#
#  Run this with /bin/sh -c <script-file-path>
#

if test -z "${SDPC_ROOT}" ; then
   echo "The SDPC environment is not set"
   exit 1
fi

# number of backup archive dbfiles to keep
archive_dbfile_num_backups=7

wait_third=1200
wait_hour=3600
wait_day=86400
wait_week=604800

# Schedule daily tasks at 4am:
daily_time="4am today"

PGMNAME="$(basename $0)"

svc_up()
{
   $SDPC_S6_ROOT/bin/s6-svstat -o up $SDPC_PIPE_DIR/services/$1
}

trace_message()
{
   args="$@"
   #echo "$args"
}

expire_dir_files()
{
   mmin_arg="$1"
   shift

   dir="$1"
   shift

   other_args="$@"

   if ! test -d $dir ; then
      echo "PGMNAME: WARNING: nonexistent directory: $dir"
      return
   fi

   find_args="$dir -mindepth 1 -maxdepth 1 -type f -mmin $mmin_arg $other_args"

   paths=$(find $find_args)
   if ! test x"$paths" = x ; then
      find $find_args -delete
   fi

   echo "$PGMNAME: cleanup $dir"
}

keep_max_files()
{
   num_keep="$1"
   path="$2"

   dir=$(dirname $path)
   name=$(basename $path)

   # 1) list matching files order of modification time
   # 2) drop the last num_keep list entries
   # 3) delete the remaining files, if any
   files_to_delete=$(find $dir -mindepth 1 -maxdepth 1 -type f -name "${name}.*" -printf '%T@\t%p\n' | \
                     sort -t $'\t' -g | \
                     head -n -$num_keep | \
                     cut -d $'\t' -f 2-)

   if test x"$files_to_delete" != x ; then
      /bin/rm -f $files_to_delete
   fi
}

last_archive_backup_t=0
backup_archive_dbfile()
{
    archive_bkp_dir="$1"

    if ! test -f $SDPC_ARCHIVE_DBFILE ; then
       return
    fi

    last_mod_t=$(stat -c %Y $SDPC_ARCHIVE_DBFILE)
    if test $(($last_mod_t - $last_archive_backup_t)) -le 0 ; then
       # Do nothing if the db file hasn't changed since the last backup
       return
    fi

    echo "$PGMNAME: backing up $SDPC_ARCHIVE_DBFILE"
    last_archive_backup_t=$last_mod_t

    if ! test -d $archive_bkp_dir ; then
       mkdir -p $archive_bkp_dir
    fi

    sqlbkp="${archive_bkp_dir}/$(basename $SDPC_ARCHIVE_DBFILE)"
    suffix="$(date -u +%Y%m%dT%H%M%SZ)"
    sqlite3 -cmd ".timeout 10000" $SDPC_ARCHIVE_DBFILE ".backup ${sqlbkp}.${suffix}"

    keep_max_files $archive_dbfile_num_backups $sqlbkp
}

replace_old_subdirs_with_tarfiles()
{
  mmin_arg="$1"
  path="$2"

  dirlist=$(find $path -mindepth 2 -maxdepth 2 -type d -name "???" -mmin $mmin_arg)
  if test -z "$dirlist" ; then
     return
  fi

  echo "$PGMNAME: tar old subdirs: $path"

  for dir in $dirlist ; do
     parent_dir="$(dirname $dir)"
     subdir="$(basename $dir)"
     tar czf "${dir}.tar.gz" --remove-files -C "$parent_dir" "$subdir"
  done
}

do_third()
{
  trace_message third

  if test $SDPC_RADREF_CRON_TRIGGER_LEVEL2 -ne 0 ; then
     flush_radref_pending.sh $SDPC_RADREF_SEARCH
  fi
}

perform_inr_qa_check()
{
   d_path=$1

   # change name to avoid a race condition
   d_path_working=${d_path}.working
   /bin/mv $d_path $d_path_working

   date_string=$(basename $d_path)
   _msg=$(inr_qa_cron.sh $date_string)

   if test $? -eq 0 ; then
      /bin/rm -f $d_path_working
   else
      /bin/mv $d_path_working ${d_path}.failed
   fi

   if test -n "$_msg" ; then
      _tstamp=$(date +%Y-%m-%dT%H:%M:%S-%Z)
      echo "$_tstamp:$_msg" >> $SDPC_PIPE_DIR/inr/quality/inr_qa.log 2>&1
   fi
}
export -f perform_inr_qa_check

manage_inr_qa_checking()
{
  if test $SDPC_INR_QA_ENABLE -eq 0 ; then
     return
  fi

  root_work_dir="$SDPC_PIPE_DIR/inr/quality"
  if ! test -d $root_work_dir ; then
     return
  fi

  # Do we have any dates to process?
  date_file_paths=$(find $root_work_dir -mindepth 1 -maxdepth 1 -type f -name "????-??-??" | sort)
  if test -z "$date_file_paths" ; then
     return
  fi
  num_dates=$(echo $date_file_paths | wc -w)

  tstamp_file="$root_work_dir/timestamp"
  # The timestamp file always contains the time_t value
  # for when the INR QA script was last run.

  timer_expired=0
  if test -f $tstamp_file ; then
     timet_prev=$(cat $tstamp_file)
     timet_now=$(date +%s)
     if test $(($timet_now - $timet_prev)) -gt 86400 ; then
	timer_expired=1
     fi
  fi

  if test $num_dates -lt 2 -a $timer_expired -eq 0 ; then
      # Timer hasn't expired, and no backlog, so do nothing and return.
      return
  fi

  # Update the timestamp file before processing
  if test -n "$SDPC_INR_QA_TIME" ; then
     date --date "$SDPC_INR_QA_TIME today" +%s > $tstamp_file
  else
     date +%s > $tstamp_file
  fi

  # Process dates serially
  for d in $date_file_paths ; do
      perform_inr_qa_check $d
  done
}

do_hourly()
{
  trace_message hourly
  manage_inr_qa_checking
}

do_daily()
{
  trace_message daily

  mmin_arg="+1440"

  isup=$(svc_up register)
  if test x"$isup" = xtrue ; then
     # backup archive sqlite database file
     backup_archive_dbfile "$SDPC_ARCHIVE_DIR/backup"
  fi

  # delete INR GUI directory files
  isup=$(svc_up inr)
  if test x"$isup" = xtrue ; then
     gui_dir="$SDPC_PIPE_DIR/stage/granules/inr_output/monitor_cache"
     if test -d $gui_dir ; then
        expire_dir_files $mmin_arg $gui_dir
     fi
  fi

  # delete EMPTY slurm batch log files
  slurm_log_services=("level1a" "level1b" "level2" "level3" "level1b_nrt1" "level1b_nrt2" "level2_nrt")
  for d in ${slurm_log_services[@]}; do
      isup=$(svc_up $d)
      if test x"$isup" = xtrue ; then
         expire_dir_files $mmin_arg "$SDPC_PIPE_DIR/log/$d/slurm" "-size 0"
      fi
  done

  # pack old daily ASDC push/pull history into tar files
  isup=$(svc_up asdc)
  if test x"$isup" = xtrue ; then
     asdc_dir="$SDPC_ARCHIVE_DIR/asdc"
     if test -d $asdc_dir ; then
        replace_old_subdirs_with_tarfiles "+2880" $asdc_dir
     fi
  fi

  # Expire old data products in the NRT archive
  num_days_kept=$SDPC_NRT_NUM_DAYS_KEPT
  if test $num_days_kept -gt 0 ; then
     expire_nrt.sh $num_days_kept
  fi
}

do_weekly()
{
  trace_message weekly
  mmin_arg="+10080"
  isup=$(svc_up inr)
  if test x"$isup" = xtrue ; then
     expire_dir_files $mmin_arg "$SDPC_PIPE_DIR/inr/scantailoring"
  fi
}

Sleep_Pid=0

handle_signal()
{
  echo "$PGMNAME: caught signal: exiting"
  if test $Sleep_Pid -ne 0 ; then
     # If the sleep process still exists, kill it.
     kill -0 $Sleep_Pid && kill -9 $Sleep_Pid
  fi
  exit
}
trap handle_signal SIGINT SIGTERM

sleep_minutes()
{
   num_minutes=$1

   # To ensure signals can interrupt this without waiting
   # for sleep to complete, run sleep in the background
   # and wait for it. Bash 'wait' will catch the signal,
   # and immediately run the corresponding trap function.

   minutes=$(seq $num_minutes)
   for i in $minutes ; do
      sleep 60 &
      Sleep_Pid=$!
      wait $Sleep_Pid
   done
}

usage_message()
{
   echo "Usage: $(basename $0) [options]"
   echo "   --help     Print this usage message"
   echo "   --service  Run as a service"
   echo 0
}

main()
{
   if test $# -ne 0 ; then
      case "$1" in
         --service) service=yes
           ;;
         --help) usage_message
           ;;
         *) usage_message
           ;;
      esac
   fi

   echo "$PGMNAME: started"

   while true ; do

      now=$(date +%s)

      if test $now -gt $(($last_third + $wait_third)) ; then
         do_third
	 last_third=$now
      fi

      if test $now -gt $(($last_hourly + $wait_hour)) ; then
         do_hourly
	 last_hourly=$now
      fi

      if test $now -gt $(($last_daily + $wait_day)) ; then
         do_daily
	 last_daily=$(date --date "$daily_time" +%s)
      fi

      if test $now -gt $(($last_weekly + $wait_week)) ; then
         do_weekly
	 last_weekly=$now
      fi

      if test -z "$service" ; then
         break
      fi

      sleep_minutes 5

   done

   echo "$PGMNAME: exiting"
}

main "$@"
