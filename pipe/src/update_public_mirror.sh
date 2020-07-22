#! /bin/sh

set -u
set -e

: "${SDPC_RUN_DIR_MASTER:?SDPC_RUN_DIR_MASTER not set}"

_source_root_dir="${SDPC_RUN_DIR_MASTER}/public_mirror"

# Silently do nothing if the top-level source directory doesn't exist
if ! test -d $_source_root_dir ; then
   exit 0
fi

# If SDPC_RSYNC_TEST is defined, use the --dry-run option"
: "${SDPC_RSYNC_TEST:=OFF}"
if test x"$SDPC_RSYNC_TEST" != xOFF ; then
   rsync_test="--dry-run"
else
   rsync_test=""
fi

PUBLIC_HOST=waps.cfa.harvard.edu
PUBLIC_USER=jhouck
PUBLIC_DIR=/data/www/cgi/sao_atmos/data/tempo_sdpc/
PUBLIC_URL="${PUBLIC_USER}@${PUBLIC_HOST}:${PUBLIC_DIR}"

# find top-level subdirectories of $_source_root_dir changed within the last 24hrs
_dirs=$(find $_source_root_dir -maxdepth 1 -mindepth 1 -type d -mtime -1)

# rync options:
#  -r = recursive
#  -L = follow links
#  -p = preserve permissions
#  -t = preserve modification times
#  -v = verbose
#  -e = specific ssh command to use

for _d in $_dirs; do
    rsync $rsync_test -rLptv --out-format="%l %f" -e 'ssh -o ForwardX11=no' $_d $PUBLIC_URL
done

INDEX_GEN="/data/www/cgi/sao_atmos/index_generator.py"

ssh -o ForwardX11=no ${PUBLIC_USER}@${PUBLIC_HOST} python $INDEX_GEN $PUBLIC_DIR
