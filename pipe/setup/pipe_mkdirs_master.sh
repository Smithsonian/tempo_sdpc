#! /bin/sh

set -u
set -e

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"

if test -e $SDPC_RUN_DIR_MASTER ; then
   echo "File exists: $SDPC_RUN_DIR_MASTER"
   exit 1
fi

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set, run this with sdpcrun.sh}"

if ! test -e $SDPC_ANCILLARY_ROOT ; then
   echo "Directory does not exist: $SDPC_ANCILLARY_ROOT"
   exit 1
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

PIPE_DIRS="incoming/l0ccd incoming/telem \
           L0/repro L1/repro \
           L2/repro L2/incoming L2/inputs L2/entry \
	   log"

mkdirlist $SDPC_RUN_DIR_MASTER "$PIPE_DIRS"

ln -s -n $SDPC_ANCILLARY_ROOT $SDPC_RUN_DIR_MASTER/ancillary
