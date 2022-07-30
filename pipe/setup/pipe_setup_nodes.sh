#! /bin/sh

set -u
set -e

: "${SDPC_PIPE_NAME:?SDPC_PIPE_NAME not set -- source sdpc_env.sh}"
: "${SDPC_ROOT:?SDPC_ROOT not set -- source sdpc_env.sh}"
: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set -- source sdpc_env.sh}"

if test $# -ne 1 ; then
   echo "Usage:  $0  <slurm-partition-name>"
   exit 0
fi

part_name=$1

# list of compute nodes in this partition
_node_list=$(sinfo -p "$part_name" -h -N -o %n)

for node in $_node_list; do
  printf "Creating pipeline node directory: $node:$SDPC_RUN_DIR\n"
  ssh $node $SDPC_ROOT/bin/pipe_mkdirs_node.sh $SDPC_PIPE_DIR $SDPC_ROOT
done
