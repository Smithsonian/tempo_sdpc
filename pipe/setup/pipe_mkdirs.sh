#! /bin/sh

set -u
set -e

: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set, run this with sdpcrun.sh}"

if test -e $SDPC_PIPE_DIR ; then
   echo "File exists: $SDPC_PIPE_DIR"
   exit 1
fi

printf "Creating pipeline directory: $SDPC_PIPE_DIR\n"

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
           stage/granules/inr_output/monitor_cache \
           stage/granules/level2_input \
           stage/scans \
           repro/L0 repro/L1 repro/L2 \
           ctrl \
	   log"

mkdirlist $SDPC_PIPE_DIR "$PIPE_DIRS"
