#! /bin/sh

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"

set -u
set -e

_telem_dir="/data/tempo/sdpc/test_data/level0/2018jun29_fixup/incoming_newnames/telem"

assert_dir_exists()
{
   dir=$1

   if ! test -e $dir ; then
      echo "Directory does not exist: $dir"
      exit 1
   fi
}

_target_dir="$SDPC_RUN_DIR_MASTER/L0/incoming/telem"

assert_dir_exists $_telem_dir
assert_dir_exists $_target_dir

ln -s -t $_target_dir $_telem_dir/*.nc
