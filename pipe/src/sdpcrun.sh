#! /bin/sh

exit_usage()
{
   echo "Usage: $(basename $0) [options] SHELL"
   echo "   Options:"
   echo "   -i SDPC_PIPE_NAME   Set pipeline name"
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
      echo "ERROR: missing required argument"
      exit_usage 1
   fi

   env_script="$top/etc/sdpc_env.sh"
   if ! test -f $env_script ; then
      echo "ERROR:  file not found: $env_script"
      exit 1
   fi
   . $env_script

   export LD_LIBRARY_PATH="$top/lib"
   exec "$@"
}

main "$@"
