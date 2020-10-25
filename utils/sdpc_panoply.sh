#! /bin/sh

top=$(dirname $0)/..
. $top/etc/sdpc_env.sh

export JAVA_HOME="$SDPC_OTS_ROOT/jdk1.8.0_91"
export PATH="$JAVA_HOME/bin:$PATH"

if ! command -v panoply &> /dev/null ; then
   echo "*** Error: panoply not found"
   exit 1
fi

panoply "$@"
