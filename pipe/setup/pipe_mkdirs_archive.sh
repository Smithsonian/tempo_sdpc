#! /bin/sh

set -u
set -e

: "${SDPC_ARCHIVE_DIR:?SDPC_ARCHIVE_DIR not set, run this with sdpcrun.sh}"

if test -d $SDPC_ARCHIVE_DIR ; then
   printf "Archive directory exists: $SDPC_ARCHIVE_DIR\n"
else
   printf "Creating archive directory: $SDPC_ARCHIVE_DIR\n"
fi

mkdirlist()
{
   rootdir="$1"
   subdirs="$2"

   mkdir -p $rootdir
   for sub in $subdirs ; do
      mkdir -p $rootdir/$sub
   done
}

ARCHIVE_DIRS="L0 L1 L2 L3"
ARCHIVE_DIRS_NRT="NRT/L1 NRT/L2 NRT/L3"

mkdirlist $SDPC_ARCHIVE_DIR "$ARCHIVE_DIRS $ARCHIVE_DIRS_NRT"

dbfile_dir=$(dirname $SDPC_ARCHIVE_DBFILE)
if ! test -d $dbfile_dir ; then
   mkdir -p $dbfile_dir
fi
