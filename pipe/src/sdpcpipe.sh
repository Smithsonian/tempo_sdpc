#! /bin/sh
top=$(dirname $0)/..
. $top/etc/sdpc_env.sh

export PATH="$SDPC_S6_ROOT/bin:$PATH"
exec "$@"
