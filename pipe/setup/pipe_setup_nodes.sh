#! /bin/sh

set -u
set -e

if test $# -lt 3 ; then
   echo "Usage:  $0  SDPC_RUN_DIR_MASTER SDPC_ROOT [node-list]"
   exit 0
fi

_run_dir_master=$1
shift

_root_dir=$1
shift

_node_list="$@"

_pipe_name=$(basename $_run_dir_master)

export SDPC_PIPE_NAME=$_pipe_name
. $_root_dir/etc/sdpc_env.sh

for node in $_node_list; do
  printf "Creating pipeline node directory: $node:$SDPC_RUN_DIR\n"
  ssh $node $_root_dir/bin/pipe_mkdirs_node.sh $_run_dir_master $_root_dir
done
