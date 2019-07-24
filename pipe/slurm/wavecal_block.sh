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
    ;;

 * )
    adjust=""
    ;;
esac

ARGS="$adjust -g $BAND_NAME -b $block_spec -c ${etc_dir}/wavecal.cfg"

result_file $INPUT_FILE

#wavecal_driver -v $ARGS $INPUT_FILE $RESULT_FILE > ${RESULT_FILE}.log 2>&1
wavecal_driver $ARGS $INPUT_FILE $RESULT_FILE

