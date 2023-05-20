#! /bin/sh

: "${SDPC_ROOT:?SDPC_ROOT not set}"
: "${SDPC_PIPE_NAME:?SDPC_PIPE_NAME not set}"

if ! test -d "$SDPC_PIPE_DIR" ; then
   printf "*** nonexistent pipeline directory: SDPC_PIPE_DIR=$SDPC_PIPE_DIR\n"
   exit 1
fi

pipe_owner="temposdpc"
when="yesterday"
when_end="now"
date_fmt="+%Y-%m-%dT%H:%M:%S"

timet_beg=""
timet_end=""
utc_beg=""
utc_end=""
loc_beg=""
loc_end=""

make_times()
{
  when="$1"

  timet_beg=$(date --date "$when" +%s)
  timet_end=$(date --date "$when_end" +%s)

  utc_beg=$(date -u --date @$timet_beg "${date_fmt}Z")
  utc_end=$(date -u --date @$timet_end "${date_fmt}Z")

  loc_beg=$(date    --date @$timet_beg "${date_fmt}")
  loc_end=$(date    --date @$timet_end "${date_fmt}")
}

separator()
{
  printf "\n"
}

service_check()
{
  echo "Errors/warnings in service log files updated since $utc_beg:"
  logdirs=$(find $SDPC_PIPE_DIR/log -mindepth 1 -maxdepth 1 -type d)
  for d in $logdirs ; do
      strings=$(filter_s6_log --beg $timet_beg --end $timet_end -i --regex "error|warn" $d | s6-tai64nlocal)
      if test -n "$strings" ; then
         printf "*** CHECK THESE SERVICE LOG MESSAGES ***\n"
	 printf "$d:\n$strings\n"
      else
         printf "\tNONE:\t$d\n"
      fi
  done

  printf "\nRepro directory contents:\n"
  dirs=$(find $SDPC_PIPE_DIR/repro -mindepth 1 -type d)
  for d in $dirs ; do
      files=$(find $d -type f)
      if test -n "$files" ; then
         printf "*** GRANULES NEED REPROCESSING ***\n"
         echo "$files"
      else
         printf "\tNONE:\t$d\n"
      fi
  done
}

cluster_check()
{
  printf "Slurm job summary:\n"
  num_completed=$(sacct -u $pipe_owner -n -A $SDPC_PIPE_NAME -S $loc_beg -E $loc_end --state CD | wc -l)
  printf "\t$num_completed:\tjobs completed\n"

  sacct_fmt="Start,Elapsed,NodeList,TotalCPU,AllocCPUS,State,ExitCode,JobID,AveVMSize,JobName%16,Comment%45"
  check_states="BF,CA,DL,F,NF,OOM,PR,RQ,RS,RV,S,TO"
  sacct_args="-u $pipe_owner -A $SDPC_PIPE_NAME -S $loc_beg -E $loc_end --format $sacct_fmt --state=$check_states"
  issues=$(sacct -n $sacct_args)

  if test -n "$issues" ; then
     num_issues=$(echo "$issues" | wc -l)
  else
     num_issues=0
  fi
  printf "\t$num_issues:\tjobs with unexpected states\n"

  nonempty_batch_logs=$(find $SDPC_PIPE_DIR/log -type f -path "*/slurm/*" -size +0c)
  if test -n "$nonempty_batch_logs" ; then
     printf "*** CHECK THESE LOGS ***\n"
     echo "$nonempty_batch_logs"
  else
     printf "\tOK:\tbatch logs\n"
  fi

  if test -n "$issues" ; then
     printf "unexpected job states: "
     printf "*** CHECK THESE JOBS ***\n"
     issues_with_header=$(sacct $sacct_args)
     echo "$issues_with_header"
  fi
}

show_status()
{
  echo "#   SDPC_PIPE_NAME:  $SDPC_PIPE_NAME"
  echo "# start time: $utc_beg"
  echo "#   end time: $utc_end"
  cluster_check
  separator
  service_check
  separator
  asdc_upload_report.py --start $utc_beg --end $utc_end
}

exit_usage()
{
  status=$1
  printf "Usage: %s [options]\n" $(basename $0)
  printf "Options:\n"
  printf "   --help          Print this usage message\n"
  printf "   --user USER     Target pipeline owner (default=temposdpc)\n"
  printf "   --start WHEN    Start time to be interpreted by 'date'\n"
  printf "   --end WHEN      End time to be interpreted by 'date'\n"
  exit $status
}

main()
{
   # Process optional args
   while [ "$#" != "0" ]
   do
     case "$1" in
       --*)
         case "$1" in
           --help)
             exit_usage 0
             ;;
           --user)
              shift;
              pipe_owner="$1"
              shift;
              ;;
           --start)
              shift;
              when="$1"
              shift;
              ;;
           --end)
              shift;
              when_end="$1"
              shift;
              ;;
           *)
             echo "Unknown option: $1"
             exit 1
             ;;
         esac
         ;;
       *)
         break
         ;;
     esac
   done

   if test "$#" -ne 0 ; then
      echo "*** Error: unexpected command line arguments"
      echo "$@"
      exit_usage 1
   fi

   make_times "$when"
   show_status
}

main "$@"
