#! /bin/sh

set -u
set -e

: "${SDPC_PIPE_NAME:?SDPC_PIPE_NAME not set -- source sdpc_env.sh}"
: "${SDPC_ROOT:?SDPC_ROOT not set -- source sdpc_env.sh}"
: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set -- source sdpc_env.sh}"

exit_usage()
{
   echo "Usage: $(basename $0) [options]"
   echo "   Options:"
   echo "   --help     Print this listing"
   echo "   --config   Pipeline context: live | cache | repro0 | repro1"
   exit "$1"
}

list_up=""
list_down=""

service_states_for_context()
{
   context="$1"
   ioc_srvs="level0"
   lev1a_srvs="level1a inr trend"
   always_up="level1b level2 level3 register pipecron"
   always_down="asdc iocpull iocpullraw"
   case "$context" in
      live)
      list_up="$ioc_srvs $lev1a_srvs $always_up"
      list_down="$always_down"
      ;;
      cache)
      list_up="$lev1a_srvs $always_up"
      list_down="$always_down $ioc_srvs"
      ;;
      repro0)
      list_up="$lev1a_srvs $always_up"
      list_down="$always_down $ioc_srvs"
      ;;
      repro1)
      list_up="$always_up"
      list_down="$always_down $ioc_srvs $lev1a_srvs"
      ;;
      *)
      echo "*** Error: unsupported pipeline context: $context"
      exit 1
      ;;
   esac
}

service_default_up()
{
   list="$1"
   if test x"$list" = x ; then
      return
   fi
   for srv in $list ; do
       /bin/rm -f $SDPC_PIPE_DIR/services/$srv/down
   done
   echo "Services up: $list"
}

service_default_down()
{
   list="$1"
   if test x"$list" = x ; then
      return
   fi
   for srv in $list ; do
       touch $SDPC_PIPE_DIR/services/$srv/down
   done
   echo "Services down: $list"
}

# Process optional args
while [ "$#" != "0" ]
do
   case "$1" in
     --*)
       case "$1" in
         --help)
           exit_usage 0
           ;;
         --config)
           shift
           if test $# = 0 ; then
              echo "*** Error: missing config argument"
              exit 1
           fi
           service_states_for_context "$1"
           shift
           ;;
         --*)
           exit_usage 0
           ;;
       esac
       ;;
     *)
       exit_usage 0
     ;;
   esac
done

pipe_mkdirs_archive.sh
pipe_mkdirs_master.sh

inr_mkdirs.sh
inr_config.sh

/bin/cp -r $SDPC_ROOT/etc $SDPC_PIPE_DIR
/bin/mv $SDPC_PIPE_DIR/etc/services $SDPC_PIPE_DIR

service_default_up "$list_up"
service_default_down "$list_down"

if test -f $SDPC_ARCHIVE_DBFILE ; then
   echo "WARNING: file exists: $SDPC_ARCHIVE_DBFILE"
fi
