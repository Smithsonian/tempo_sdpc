#
#  Run this with /bin/sh -c <script-file-path>
#

if test -z "${SDPC_ROOT}" ; then
   echo "The SDPC environment is not set"
   exit 1
fi

# number of backup archive dbfiles to keep
archive_dbfile_num_backups=7

wait_hour=3600
wait_day=86400
wait_week=604800

PGMNAME="$(basename $0)"

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
    sqlite3 $SDPC_ARCHIVE_DBFILE ".backup ${sqlbkp}.${suffix}"

    keep_max_files $archive_dbfile_num_backups $sqlbkp
}

do_hourly()
{
  trace_message hourly
}

do_daily()
{
  trace_message daily

  mmin_arg="+1440"

  # backup archive sqlite database file
  backup_archive_dbfile "$SDPC_ARCHIVE_DIR/backup"

  # delete INR GUI directory files
  gui_dir="$SDPC_RUN_DIR_MASTER/stage/granules/inr_output/GUI"
  if test -d $gui_dir ; then
     expire_dir_files $mmin_arg $gui_dir
  fi

  # delete EMPTY slurm batch log files
  slurm_log_dirs=("level1a" "level1b" "level2")
  for d in ${slurm_log_dirs[@]}; do
      expire_dir_files $mmin_arg "$SDPC_RUN_DIR_MASTER/log/$d/slurm" "-size 0"
  done
}

do_weekly()
{
  trace_message weekly
  mmin_arg="+10080"
  expire_dir_files $mmin_arg "$SDPC_RUN_DIR_MASTER/inr/scantailoring"
}

handle_signal()
{
  echo "$PGMNAME: caught signal: exiting"
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
      wait
   done
}

main()
{
   echo "$PGMNAME: started"

   while true ; do

      now=$(date +%s)

      if test $now -gt $(($last_hourly + $wait_hour)) ; then
         do_hourly
	 last_hourly=$now
      fi

      if test $now -gt $(($last_daily + $wait_day)) ; then
         do_daily
	 last_daily=$now
      fi

      if test $now -gt $(($last_weekly + $wait_week)) ; then
         do_weekly
	 last_weekly=$now
      fi

      sleep_minutes 5

   done

   echo "$PGMNAME: exiting"
}

main
