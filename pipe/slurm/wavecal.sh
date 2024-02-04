#! /bin/bash

if test $# -ne 2 ; then
   echo "Usage: wavecal.sh FILE NUM"
   exit 1
fi

input_file=$1
num=$2

bn=$(basename $input_file .nc)

result_dir="wavecal_${bn}"
mkdir -p $result_dir  || exit 1

wrap_wavecal()
{
   band_name=$1
   block_spec=$2

   result_file="${result_dir}/${bn}_wavecal_${band_name}_${block_spec}.nc"

   etc_dir="$SDPC_PIPE_DIR/etc"

   case "$bn" in
    *RAD* )
       #adjust="--adjust"
       adjust=""
       block_args="--s_block $block_spec"
       config="${etc_dir}/wavecal_rad.cfg"
       ;;

    *IRR* )
       adjust=""
       block_args="--x_block $block_spec"
       config="${etc_dir}/wavecal_irr.cfg"
       ;;
   esac

   args="$adjust $block_args -g $band_name -c $config"

   srun --nodes 1 --ntasks 1 --exclusive \
        --job-name="wavecal" \
        --output=${result_file}.log \
         wavecal_driver $args $input_file $result_file
}
export -f wrap_wavecal
export input_file
export bn
export result_dir

parallel --delay .2 -j $SLURM_NTASKS --joblog wavecal_joblog.out \
       wrap_wavecal {1} {2}:$num ::: band_290_490_nm band_540_740_nm ::: $(seq 0 $((num-1)))

### Uncomment these two lines to collect stderr/stdout logs from wavecal_driver
wavecal_log_files=$(find $result_dir -mindepth 1 -maxdepth 1 -name "*.log")
if test x"$wavecal_log_files" != x ; then
   tar cz --remove-files -f wavecal_logs.tar.gz $wavecal_log_files
fi

case "$bn" in
 *IRR* )
    args="--meta"
    ;;

 * )
    args=""
    ;;
esac

srun --nodes 1 --ntasks 1 --output=log_wavecal_merge.txt \
     --job-name="wavecal_merge"  \
     wavecal_merge $args --delete -t $input_file $result_dir
