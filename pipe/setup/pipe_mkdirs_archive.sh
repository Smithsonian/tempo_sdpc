#! /bin/sh

set -u
set -e

: "${SDPC_ARCHIVE_DIR:?SDPC_ARCHIVE_DIR not set, run this with sdpcrun.sh}"

if test -e $SDPC_ARCHIVE_DIR ; then
   echo "File exists: $SDPC_ARCHIVE_DIR"
   exit 1
fi

mkdirlist()
{
   rootdir="$1"
   subdirs="$2"

   echo mkdir -p $rootdir
   for sub in $subdirs ; do
      echo mkdir -p $rootdir/$sub
   done
}

ARCHIVE_DIRS="L0/incoming L1/incoming L2/incoming L3/incoming"

mkdirlist $SDPC_ARCHIVE_DIR "$ARCHIVE_DIRS"
