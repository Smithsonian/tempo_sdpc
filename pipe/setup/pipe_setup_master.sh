#! /bin/sh

set -u
set -e

: "${SDPC_PIPE_NAME:?SDPC_PIPE_NAME not set -- source sdpc_env.sh}"
: "${SDPC_ROOT:?SDPC_ROOT not set -- source sdpc_env.sh}"
: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set -- source sdpc_env.sh}"

pipe_mkdirs_master.sh
pipe_mkdirs_archive.sh

inr_mkdirs.sh
inr_config.sh

/bin/cp -r $SDPC_ROOT/etc $SDPC_RUN_DIR_MASTER
/bin/mv $SDPC_RUN_DIR_MASTER/etc/services $SDPC_RUN_DIR_MASTER

