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
   echo "   --help          Print this listing"
   echo "   --config        Pipeline context: live | cache | repro0 | repro1"
   echo "   --archive NAME   Archiving will be handled by pipeline NAME"
   exit "$1"
}

list_up=""
list_down=""

service_states_for_context()
{
   context="$1"
   ioc_srvs="level0"
   lev1a_srvs="level1a inr trend"
   nrt_srvs="level1b_nrt1 level1b_nrt2 level2_nrt"
   always_up="level1b level2 level3 register pipecron"
   always_down="asdc iocpull iocpullraw"
   case "$context" in
      live)
      list_up="$ioc_srvs $lev1a_srvs $nrt_srvs $always_up"
      list_down="$always_down"
      ;;
      cache)
      list_up="$lev1a_srvs $nrt_srvs $always_up"
      list_down="$always_down $ioc_srvs"
      ;;
      repro0)
      list_up="$lev1a_srvs $always_up"
      list_down="$always_down $ioc_srvs $nrt_srvs"
      ;;
      repro1)
      list_up="$always_up"
      list_down="$always_down $ioc_srvs $lev1a_srvs $nrt_srvs"
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

make_archive()
{
   pipe_mkdirs_archive.sh
   if test -f $SDPC_ARCHIVE_DBFILE ; then
      echo "WARNING: file exists: $SDPC_ARCHIVE_DBFILE"
   fi
}

setup_archiver()
{
   archiver_name="$1"

   echo "Using external archiver: $archiver_name"
   echo "Updating service defaults:"
   service_default_down "level3 trend register"

   # We can create the archiver script to setup future shells, but
   # environment variables in the current shell must be set manually:
   _archive_dir="$(echo $SDPC_ARCHIVE_DIR | sed -e s,/$SDPC_PIPE_NAME,/$archiver_name,)"
   _dbfile="$(echo $SDPC_ARCHIVE_DBFILE | sed -e s,/$SDPC_PIPE_NAME,/$archiver_name,)"
   archiver_script="$SDPC_PIPE_DIR/etc/archiver.sh"
   echo "export SDPC_ARCHIVE_DIR=$_archive_dir" > $archiver_script
   echo "export SDPC_ARCHIVE_DBFILE=$_dbfile" >> $archiver_script
   echo "Saved archive paths in: $archiver_script"

   printf "\n***IMPORTANT: Please set these environment variables now:\n\n"
   cat $archiver_script
}

if test $# -eq 0 ; then
   exit_usage 0
fi

context=""
archiver_name=""
create_archive=yes

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
           context="$1"
           service_states_for_context "$1"
           shift
           ;;
         --archive)
           shift
           if test $# = 0 ; then
              echo "*** Error: missing config argument"
              exit 1
           fi
           archiver_name="$1"
           create_archive=no
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

pipe_mkdirs.sh

inr_mkdirs.sh
inr_config.sh

/bin/cp -r $SDPC_ROOT/etc $SDPC_PIPE_DIR
tar -C $SDPC_PIPE_DIR -xf $SDPC_ROOT/services.tar

echo "Configuring services for $context context:"
service_default_up "$list_up"
service_default_down "$list_down"

if test x"$create_archive" = xyes ; then
   make_archive
elif test -n "$archiver_name" ; then
   setup_archiver "$archiver_name"
fi
