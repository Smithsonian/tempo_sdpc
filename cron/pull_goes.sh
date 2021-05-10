#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

tmpfile=$(mktemp -p $SDPC_ANCILLARY_ROOT/src)

cleanup()
{
   /bin/rm -f $tmpfile
}
trap cleanup EXIT

rootdir="${SDPC_ANCILLARY_ROOT}/goes"
subdir="$(date +%Y/%j)"

target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

sed -e s,@TARGET_DIR@,$target_dir,g \
    $SDPC_ANCILLARY_ROOT/src/pda_lftp.script > $tmpfile

lftp -f $tmpfile

# The 'today' symlink always points to today's GOES data
cd $rootdir
ln -nf -s $subdir today
