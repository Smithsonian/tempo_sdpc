#! /bin/sh

set -u
set -e

if test $# -lt 3 ; then
   echo "Usage:  $0  SDPC_PIPE_NAME SDPC_ROOT [node-list]"
   exit 0
fi

_pipe_name=$1
shift

_root_dir=$1
shift

_node_list="$@"

export SDPC_PIPE_NAME=$_pipe_name
. $_root_dir/etc/sdpc_env.sh

for node in $_node_list; do
  ssh $node $_root_dir/bin/pipe_mkdirs_node.sh $_pipe_name $_root_dir
done
