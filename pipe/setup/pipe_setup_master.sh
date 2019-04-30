#! /bin/sh

set -u
set -e

if test $# -ne 2 ; then
   echo "Usage:  $0 SDPC_PIPE_NAME SDPC_ROOT"
   exit 0
fi

_pipe_name=$1
_root_dir=$2

export SDPC_PIPE_NAME=$_pipe_name

. $_root_dir/etc/sdpc_env.sh

pipe_mkdirs_master.sh
pipe_mkdirs_archive.sh

inr_mkdirs.sh
inr_config.sh
inr_refdata.sh

test_telem.sh

# FIXME - During operations, cron jobs will do this.
#         For testing purposes, do it here.
filedb -c $SDPC_ROOT/etc/filedb.cfg met --update
filedb -c $SDPC_ROOT/etc/filedb.cfg snow --update

/bin/cp -r $SDPC_ROOT/etc/services $SDPC_RUN_DIR_MASTER

