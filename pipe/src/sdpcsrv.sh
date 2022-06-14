#!/bin/sh

exit_usage()
{
   echo "Usage: $(basename $0) [options] {up|down|stop|cont|usr1} SERVICE"
   echo "   Options:"
   echo "   --help              Print this listing"
   echo "   --dir SERVICE_DIR   Pipeline service directory path"
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

main ()
{
   if [ "$#" = "0" ]; then
     exit_usage 0
   fi

   if test -z "$SDPC_RUN_DIR_MASTER" ; then
      _scan_dir=""
   else
      _scan_dir="$SDPC_RUN_DIR_MASTER/services"
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
           --dir)
              shift;
              _scan_dir="$1"
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

  if test x"$_scan_dir" = x ; then
    echo "SDPC_RUN_DIR_MASTER env var is undefined. Consider using the --dir option."
    exit 1
  fi

  if test "$#" -ne 2 ; then
    exit_usage 1
  fi

  _action="$1"
  _service_dir="$_scan_dir/$2"

  check_dir "$_service_dir"

  case "$_action" in
      up)
        s6-svc -u $_service_dir
      ;;

      down)
        s6-svc -d $_service_dir
      ;;

      stop)
        s6-svc -p $_service_dir
      ;;

      cont)
        s6-svc -c $_service_dir
      ;;

      usr1)
        s6-svc -1 $_service_dir
      ;;

      *)
        echo "Unknown action: $_action"
        exit_usage 1
      ;;
  esac
}

check_config
main "$@"

