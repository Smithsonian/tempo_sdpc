#! /bin/sh
#SBATCH --output=/dev/null

input_file="$1"
result_dir="$2"

wavecal_merge --delete -t $input_file $result_dir
