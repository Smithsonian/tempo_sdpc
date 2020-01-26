#! /bin/sh

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"

set -u
set -e

#_src_dir="/data/tempo/sdpc/test_data/level0/2018jun29_fixup/incoming_newnames/telem"

_src_dir="/data/tempo/sdpc/test_data/level0/2019jul21/v2/telem"
_target_dir="$SDPC_ARCHIVE_DIR/L0/D4944"

assert_dir_exists()
{
   dir=$1

   if ! test -e $dir ; then
      echo "Directory does not exist: $dir"
      exit 1
   fi
}

assert_dir_exists $_src_dir

for d in HK IRU SMC ; do
   mkdir -p $_target_dir/$d
   ln -s -t $_target_dir/$d $_src_dir/TEMPO_${d}*.nc
done

# another stupid hack (because the DRK and IRR files have the wrong times,
# and the filenames don't match the 'time' array in the file either. Sigh).
_target_dir="$SDPC_ARCHIVE_DIR/L0/D7134/HK"
mkdir -p $_target_dir
ln -s -t $_target_dir $_src_dir/TEMPO_HK_L0_V00_20190715T051005Z.nc
ln -s -t $_target_dir $_src_dir/TEMPO_HK_L0_V00_20190715T052010Z.nc

_target_dir="$SDPC_ARCHIVE_DIR/L0/D7135/HK"
mkdir -p $_target_dir
ln -s -t $_target_dir $_src_dir/TEMPO_HK_L0_V00_20190715T081005Z.nc
