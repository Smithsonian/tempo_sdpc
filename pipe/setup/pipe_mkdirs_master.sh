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

PIPE_DIRS="L0/repro L0/incoming/ccd L0/incoming/telem \
           L1/repro \
	   L2/repro L2/incoming \
	   status"

mkdirlist $SDPC_RUN_DIR_MASTER "$PIPE_DIRS"

ln -s -n $SDPC_ANCILLARY_ROOT $SDPC_RUN_DIR_MASTER/ancillary
