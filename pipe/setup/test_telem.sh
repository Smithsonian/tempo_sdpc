#! /bin/sh

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"

set -u
set -e

_src_dir="/data/tempo/sdpc/test_data/level0/2018jun29_fixup/incoming_newnames/telem"

assert_dir_exists()
{
   dir=$1

   if ! test -e $dir ; then
      echo "Directory does not exist: $dir"
      exit 1
   fi
}

assert_dir_exists $_src_dir

_target_dir="$SDPC_ARCHIVE_DIR/L0/1/4943"

for d in hk iru smc ; do
   mkdir -p $_target_dir/$d
   ln -s -t $_target_dir/$d $_src_dir/TEMPO_${d}*.nc
done
