#! /bin/sh

set -u
set -e

if test $# -ne 2 ; then
   echo "Usage:  $0  SDPC_PIPE_DIR SDPC_ROOT"
   exit 0
fi

_pipe_dir=$1
_root_dir=$2

_pipe_name=$(basename $_pipe_dir)

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

#assert_dir_absent $SDPC_NODE_DIR
assert_dir_exists $SDPC_REFDATA_DIR
assert_dir_exists $SDPC_PIPE_DIR
assert_dir_exists $SDPC_PIPE_DIR/inr

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
  ln -s -f -n $from $to
}

PIPE_NODE_DIRS="L0 L1 L2"

mkdirlist $SDPC_NODE_DIR "$PIPE_NODE_DIRS"
