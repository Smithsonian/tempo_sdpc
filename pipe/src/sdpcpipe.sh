#! /bin/sh

if test -z "$SDPC_PIPE_NAME" ; then
   printf "*** sdpcpipe.sh: SDPC_PIPE_NAME not set\n"
   exit 1
fi

top=$(dirname $0)/..
. $top/etc/sdpc_env.sh

export PATH="$SDPC_S6_ROOT/bin:$PATH"
exec "$@"
