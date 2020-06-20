#! /bin/sh
#SBATCH --output=/dev/null

input_file="$1"
result_dir="$2"

bn=$(basename $input_file)
case "$bn" in
 *IRR* )
    args="--meta"
    ;;

 * )
    args=""
    ;;
esac

### Uncomment these two lines to collect stderr/stdout logs from wavecal_driver
wavecal_log_files=$(ls $result_dir/*.log)
tar czv --remove-files -f wavecal_logs.tar.gz $wavecal_log_files

wavecal_merge $args --delete -t $input_file $result_dir
