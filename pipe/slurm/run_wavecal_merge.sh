#! /bin/sh
#SBATCH --output=/dev/null

input_file="$1"
result_dir="$2"

bn=$(basename $input_file)
case "$bn" in
 *irr* )
    args="--meta"
    ;;

 * )
    args=""
    ;;
esac

wavecal_merge $args --delete -t $input_file $result_dir
