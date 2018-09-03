#! /bin/sh

: ${SDPC_RUN_DIR?}
path="${SDPC_RUN_DIR}/ancillary/met"
met_file=`readlink $path/latest_met`
echo "$path/$met_file"
