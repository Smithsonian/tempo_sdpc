#
#  Run this with /bin/sh -c <script-file-path>
#

if test -z "${SDPC_ROOT}" ; then
   echo "The SDPC environment is not set"
   exit 1
fi

delta_sec=60
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
   mtime_arg="$1"
   shift

   dir="$1"
   shift

   other_args="$@"

   if ! test -d $dir ; then
      echo "PGMNAME: WARNING: nonexistent directory: $dir"
      return
   fi

   find_args="$dir -mindepth 1 -maxdepth 1 -type f -mtime $mtime_arg $other_args"

   paths=$(find $find_args)
   if ! test x"$paths" = x ; then
      find $find_args -delete
   fi

   echo "$PGMNAME: cleanup $dir"
}

do_hourly()
{
  trace_message hourly
}

do_daily()
{
  trace_message daily

  mtime_arg="+1"

  # delete INR GUI directory files
  expire_dir_files $mtime_arg "$SDPC_RUN_DIR_MASTER/stage/granules/inr_output/GUI"

  # delete EMPTY slurm batch log files
  slurm_log_dirs=("level1a" "level1b" "level2")
  for d in ${slurm_log_dirs[@]}; do
      expire_dir_files $mtime_arg "$SDPC_RUN_DIR_MASTER/log/$d/slurm" "-size 0"
  done
}

do_weekly()
{
  trace_message weekly
  mtime_arg="+7"
  expire_dir_files $mtime_arg "$SDPC_RUN_DIR_MASTER/inr/scantailoring"
}

main()
{
   while true ; do

      trace_message sleeping $delta_sec
      sleep $delta_sec
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

   done
}

main
