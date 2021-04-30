#! /bin/sh
#SBATCH --output=/dev/null

# Usage:
#    sbatch --array=0-9 ./wavecal_block.sh <band-name> <input-file> <result-dir>

set -e
set -u

if ! test $# -eq 3 ; then
   echo "Usage: $0 <band-name> <input-file> <result-dir>"
   exit 0
fi

BAND_NAME="$1"
INPUT_FILE="$2"
RESULT_DIR="$3"

this_block="$SLURM_ARRAY_TASK_ID"
num_blocks="$SLURM_ARRAY_TASK_COUNT"

block_spec="${this_block}:${num_blocks}"

result_file()
{
  path=$1
  bn=$(basename $path .nc)
  RESULT_FILE="${RESULT_DIR}/${bn}_wavecal_${BAND_NAME}_${block_spec}.nc"
}

etc_dir="$SDPC_ROOT/etc"

bn=$(basename $INPUT_FILE)
case "$bn" in
 *RAD* )
    adjust="--adjust"
    block_args="--s_block $block_spec"
    config="${etc_dir}/wavecal_rad.cfg"
    ;;

 *IRR* )
    adjust=""
    block_args="--x_block $block_spec"
    config="${etc_dir}/wavecal_irr.cfg"
    ;;
esac

ARGS="$adjust $block_args -g $BAND_NAME -c $config"

result_file $INPUT_FILE

wavecal_driver $ARGS $INPUT_FILE $RESULT_FILE > ${RESULT_FILE}.log 2>&1

