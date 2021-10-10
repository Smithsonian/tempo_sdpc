#! /bin/sh

set -u
set -e

if test $# -ne 2 ; then
   echo "Usage:  $0  SDPC_RUN_DIR_MASTER SDPC_ROOT"
   exit 0
fi

_run_dir_master=$1
_root_dir=$2

_pipe_name=$(basename $_run_dir_master)

if test x"$_pipe_name" != x"$_run_dir_master" ; then
   _pipe_home=$(dirname $_run_dir_master)
   export SDPC_PIPE_HOME=$_pipe_home
fi

export SDPC_PIPE_NAME=$_pipe_name
. $_root_dir/etc/sdpc_env.sh

assert_dir_exists()
{
   dir=$1

   if ! test -e $dir ; then
      echo "Directory does not exist: $dir"
      exit 1
   fi
}

assert_dir_absent()
{
   dir=$1

   if test -e $dir ; then
      echo "File exists: $dir"
      exit 1
   fi
}

assert_dir_absent $SDPC_RUN_DIR
assert_dir_exists $SDPC_REFDATA_DIR
assert_dir_exists $SDPC_RUN_DIR_MASTER
assert_dir_exists $SDPC_RUN_DIR_INR

mkdirlist()
{
   rootdir="$1"
   subdirs="$2"

   mkdir -p $rootdir
   for sub in $subdirs ; do
      mkdir -p $rootdir/$sub
   done
}

link_existing_dir()
{
  from=$1
  to=$2

  assert_dir_exists $from
  ln -s -n $from $to
}

PIPE_NODE_DIRS="L0 L1 L2"

mkdirlist $SDPC_RUN_DIR "$PIPE_NODE_DIRS"

link_existing_dir $SDPC_REFDATA_DIR $SDPC_RUN_DIR/refdata
