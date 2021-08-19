#! /bin/sh

set -u
set -e

if test $# -lt 3 ; then
   echo "Usage:  $0  SDPC_RUN_DIR_MASTER SDPC_ROOT [node-list]"
   exit 0
fi

_run_dir_master=$1
shift

if ! test -d $_run_dir_master ; then
    printf "*** Cannot find directory SDPC_RUN_DIR_MASTER = $_run_dir_master\n"
    exit 1
fi

_root_dir=$1
shift

if ! test -d $_root_dir ; then
    printf "*** Cannot find directory SDPC_ROOT = $_root_dir\n"
    exit 1
fi

_node_list="$@"

_pipe_name=$(basename $_run_dir_master)
export SDPC_PIPE_NAME=$_pipe_name

env_path="$_root_dir/etc/sdpc_env.sh"
if ! test -f "$env_path" ; then
   printf "Cannot find script file: $env_path\n"
   exit 1
fi

. "$env_path"

for node in $_node_list; do
  printf "Creating pipeline node directory: $node:$SDPC_RUN_DIR\n"
  ssh $node $_root_dir/bin/pipe_mkdirs_node.sh $_run_dir_master $_root_dir
done
