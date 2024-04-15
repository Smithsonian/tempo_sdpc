#! /bin/sh

set -u
set -e

: "${SDPC_PIPE_NAME:?SDPC_PIPE_NAME not set -- source sdpc_env.sh}"
: "${SDPC_ROOT:?SDPC_ROOT not set -- source sdpc_env.sh}"
: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set -- source sdpc_env.sh}"

exit_usage()
{
  echo "Usage:  $(basename $0) [options] <slurm-partition-name>"
  echo "Optional arguments:"
  echo " --help    Print this listing"
  echo " --check   Check that node directories exist"
  exit "$1"
}

if test $# -eq 0 ; then
   exit_usage 0
fi

_check=""

# Process optional args
while [ "$#" != "0" ]
do
  case "$1" in
    --*)
      case "$1" in
        --help)
          exit_usage 0
          ;;
        --check)
          _check="$1"
          shift;
          ;;
        *)
          echo "Unknown option: $1"
          exit 1
          ;;
      esac
      ;;
    *)
      break
      ;;
  esac
done

if test $# -gt 1 ; then
   exit_usage 1
elif test $# -eq 1 ; then
   part_name=$1
else
   # use default partition name
   if test -n "$SLURM_PARTITION" ; then
      part_name="$SLURM_PARTITION"
   elif test -n "$SBATCH_PARTITION" ; then
      part_name="$SBATCH_PARTITION"
   else
      part_name=$(sinfo -o "%P" | grep '*' | tr -d '*')
   fi
fi

# list of responding compute nodes in this partition
_node_list=$(sinfo -p "$part_name" -r -h -N -o %n)

for node in $_node_list; do
  if test -z "$_check" ; then
     printf "Creating pipeline node directory: $node:$SDPC_NODE_DIR\n"
  else
     printf "Checking pipeline node directory: $node:$SDPC_NODE_DIR\n"
  fi
  ssh $node $SDPC_ROOT/bin/pipe_mkdirs_node.sh $_check $SDPC_PIPE_DIR $SDPC_ROOT
done
