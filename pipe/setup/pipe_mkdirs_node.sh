#! /bin/sh

set -u
set -e

exit_usage()
{
  echo "Usage:  $(basename $0) [options] SDPC_PIPE_DIR SDPC_ROOT"
  echo "Optional arguments:"
  echo " --help    Print this listing"
  echo " --check   Check that node directories exist"
  exit "$1"
}

if test $# -eq 0 ; then
   exit_usage 0
fi

_check=""

# Process optional args
while [ "$#" != "0" ]
do
  case "$1" in
    --*)
      case "$1" in
        --help)
          exit_usage 0
          ;;
        --check)
          _check="$1"
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

if test $# -ne 2 ; then
   exit_usage 1
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

PIPE_NODE_DIRS="L0 L1 L2"

if test -z "$_check" ; then
   mkdirlist $SDPC_NODE_DIR "$PIPE_NODE_DIRS"
else
   for sub in $PIPE_NODE_DIRS ; do
       assert_dir_exists $SDPC_NODE_DIR/$sub
   done
fi
