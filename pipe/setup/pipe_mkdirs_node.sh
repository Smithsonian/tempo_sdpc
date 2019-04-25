#! /bin/sh

set -u
set -e

if test $# -ne 2 ; then
   echo "Usage:  $0  SDPC_PIPE_NAME SDPC_ROOT"
   exit 0
fi

_pipe_name=$1
_root_dir=$2

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
assert_dir_exists $SDPC_INR_RUN_DIR

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

PIPE_NODE_DIRS="L0/incoming L1 L2"

mkdirlist $SDPC_RUN_DIR "$PIPE_NODE_DIRS"

link_existing_dir $SDPC_REFDATA_DIR $SDPC_RUN_DIR/refdata
link_existing_dir $SDPC_RUN_DIR_MASTER/ancillary $SDPC_RUN_DIR/ancillary
link_existing_dir $SDPC_RUN_DIR_MASTER/L2/incoming $SDPC_RUN_DIR/L2/incoming
link_existing_dir $SDPC_RUN_DIR_MASTER/L2/repro $SDPC_RUN_DIR/L2/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L1/repro $SDPC_RUN_DIR/L1/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L0/repro $SDPC_RUN_DIR/L0/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L0/incoming/telem $SDPC_RUN_DIR/L0/incoming/telem
link_existing_dir $SDPC_INR_RUN_DIR/Staging/Granules $SDPC_RUN_DIR/L1/radiance_inr_staging
