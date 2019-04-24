#! /bin/sh

set -u
set -e

if test $# -ne 1 ; then
   echo "Usage:  $0 <pipe_label>"
   exit 0
fi

PIPE_LABEL=$1

# Local root processing directory on each compute node:
NODE_ROOT_DIR=/scratch/$PIPE_LABEL/sdpc_run_dir

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

assert_dir_absent $NODE_ROOT_DIR

: "${SDPC_REFDATA_DIR:?SDPC_REFDATA_DIR not set, run this with sdpcrun.sh}"
assert_dir_exists $SDPC_REFDATA_DIR

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"
assert_dir_exists $SDPC_RUN_DIR_MASTER

: "${SDPC_INR_RUN_DIR:?SDPC_INR_RUN_DIR not set, run this with sdpcrun.sh}"
assert_dir_exists $SDPC_INR_RUN_DIR

mkdirlist()
{
   rootdir="$1"
   subdirs="$2"

   echo mkdir -p $rootdir
   for sub in $subdirs ; do
      echo mkdir -p $rootdir/$sub
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

mkdirlist $NODE_ROOT_DIR "$PIPE_NODE_DIRS"

link_existing_dir $SDPC_REFDATA_DIR $NODE_ROOT_DIR/refdata
link_existing_dir $SDPC_RUN_DIR_MASTER/ancillary $NODE_ROOT_DIR/ancillary
link_existing_dir $SDPC_RUN_DIR_MASTER/L2/incoming $NODE_ROOT_DIR/L2/incoming
link_existing_dir $SDPC_RUN_DIR_MASTER/L2/repro $NODE_ROOT_DIR/L2/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L1/repro $NODE_ROOT_DIR/L1/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L0/repro $NODE_ROOT_DIR/L0/repro
link_existing_dir $SDPC_RUN_DIR_MASTER/L0/incoming/telem $NODE_ROOT_DIR/L0/incoming/telem
link_existing_dir $SDPC_INR_RUN_DIR/Staging/Granules $NODE_ROOT_DIR/L1/radiance_inr_staging
