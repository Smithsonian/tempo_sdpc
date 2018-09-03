#! /bin/sh

: ${SDPC_RUN_DIR?}
path="${SDPC_RUN_DIR}/ancillary/snow"
snow_file=`readlink $path/latest_snow`
echo "$path/$snow_file"
