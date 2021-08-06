#!/bin/sh

exit_usage()
{
   echo "Usage: $(basename $0) {status|start|stop} /path/to/services"
   exit "$1"
}

check_config()
{
   if test "X${SDPC_S6_ROOT}" = "X"
   then
     echo "SDPC_S6_ROOT env var is not defined.  Run this script via sdpcrun.sh"
     exit 1
   fi
}

check_dir()
{
  if ! test -d "$1" ; then
     echo "Directory $1 does not exist"
     exit 1
  fi
}

do_start()
{
  check_dir "$1"

  echo "Starting services in $1"
  $SDPC_S6_ROOT/bin/s6-svscan "$1" &
}

do_stop()
{
  check_dir "$1"

  echo "Stopping services in $1"
  $SDPC_S6_ROOT/bin/s6-svscanctl -t "$1"
}

do_status()
{
  check_dir "$1"
  scan_dir="$1"

  service_dirs=$(find "$1" -maxdepth 1 -mindepth 1 -type d | grep -v .s6-svscan | sort)

  # Are supervisors running?
  running=""
  for srv in $service_dirs ; do
      if $SDPC_S6_ROOT/bin/s6-svok $srv ; then
         running="1"
         break
      fi
  done
  if test "X$running" == "X" ; then
     printf "Not started\n"
     exit 0
  fi

  if test "X$service_dirs" != "X" ; then
     printf "%12s \tSTATE\tLOG\n" SERVICE
  fi
  for srv in $service_dirs ; do

     # Check service state
     srv_state=$($SDPC_S6_ROOT/bin/s6-svstat -o up $srv)
     case $srv_state in
        true) srv_state=up ;;
        *) srv_state=down ;;
     esac

     log_msg=""
     if test -d $srv/log ; then
        # Check logger state
        log_state=$($SDPC_S6_ROOT/bin/s6-svstat -o up $srv/log)
        case $log_state in
           true) log_state=up ;;
           *) log_state=down ;;
        esac
        log_msg="\t${log_state}"
     fi

     printf "%12s \t${srv_state}${log_msg}\n" $(basename $srv)
  done
}

main ()
{
   if [ "$#" = "0" ]; then
     exit_usage 0
   fi

   # Process optional args
   while [ "$#" != "0" ]
   do
     case "$1" in
       --*)
         case "$1" in
           --help)
             exit_usage 0
             ;;
          # --init)
          #    _do_init=1
          #    shift;
          #    ;;
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

  if [ "$#" != "2" ]; then
    exit_usage 1
  fi

  _action="$1"
  _scan_dir="$2"

  case "$_action" in
      start)
        do_start $_scan_dir
      ;;

      stop)
        do_stop $_scan_dir
      ;;

      status)
        do_status $_scan_dir
      ;;

      *)
        echo "Unknown action: $_action"
        exit_usage 1
      ;;
  esac
}

check_config
main "$@"

