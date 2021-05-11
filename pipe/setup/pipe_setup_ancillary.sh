#! /bin/sh

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set}"

set -e
set -u

# FIXME - During operations, cron jobs will do this.
#         For testing purposes, do it here.
filedb_cfg=$SDPC_RUN_DIR_MASTER/etc/filedb.cfg
filedb -c $filedb_cfg met:synth --update
filedb -c $filedb_cfg met:hires --update
filedb -c $filedb_cfg met:lores --update
filedb -c $filedb_cfg snow --update
