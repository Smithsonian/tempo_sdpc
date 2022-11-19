#! /bin/sh

exit_usage()
{
   echo "Usage: $(basename $0) [options] command"
   echo "   Options:"
   echo "   -i SDPC_PIPE_NAME   Set pipeline name"
   echo "   -u SDPC_USER        Set pipeline owner (only for directory paths)"
   echo "   --help              Print this listing"
   exit "$1"
}

main ()
{
   top=$(dirname $0)/..

   # Process optional args
   while [ "$#" != "0" ]
   do
     case "$1" in
       -i)
         shift
         if test -z "$1" ; then
            exit_usage 1
         fi
         export SDPC_PIPE_NAME="$1"
         shift
         ;;
       -u)
         shift
         if test -z "$1" ; then
            exit_usage 1
         fi
         # don't export this
         SDPC_USER="$1"
         shift
         ;;
       --*)
         case "$1" in
           --help)
             exit_usage 0
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

   other_args="$@"
   if test -z "$other_args" ; then
      exit_usage 0
   fi

   env_script="$top/etc/sdpc_env.sh"
   if ! test -f $env_script ; then
      echo "ERROR:  file not found: $env_script"
      exit 1
   fi
   . $env_script

   archiver_script="$SDPC_PIPE_DIR/etc/archiver.sh"
   if test -f $archiver_script ; then
      echo "Getting archive paths from: $archiver_script"
      . $archiver_script
   fi

   export LD_LIBRARY_PATH="$top/lib"
   exec "$@"
}

main "$@"
