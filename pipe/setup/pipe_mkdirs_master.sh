#! /bin/sh

set -u
set -e

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set, run this with sdpcrun.sh}"

if test -e $SDPC_RUN_DIR_MASTER ; then
   echo "File exists: $SDPC_RUN_DIR_MASTER"
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

PIPE_DIRS="stage/level0/ccd \
           stage/level0/telem \
           stage/level0/spool \
           stage/granules/inr_input \
           stage/granules/inr_output \
           stage/granules/level2_input \
           stage/scans \
           repro/L0 repro/L1 repro/L2 \
	   log public_mirror"

mkdirlist $SDPC_RUN_DIR_MASTER "$PIPE_DIRS"
