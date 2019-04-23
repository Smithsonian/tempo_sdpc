#! /bin/sh
top=$(dirname $0)/..
. $top/etc/sdpc_env.sh

export LD_LIBRARY_PATH="$top/lib"
exec "$@"
