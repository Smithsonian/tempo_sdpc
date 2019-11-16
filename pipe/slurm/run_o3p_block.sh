#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

run_dir=$1

block_run_subdir=$(printf "block_%03d" $SLURM_ARRAY_TASK_ID)
cd "${run_dir}/O3PROF/${block_run_subdir}"

# Load default config parameters
config_file="$SDPC_ROOT/etc/o3_profile/o3_profile.rc"
. $config_file

export PGSMSG="${SDPC_ROOT}/msgs"
export PGS_PC_INFO_FILE="$pcf_file"

srun --ntasks=1 --cpus-per-task=1 --exclusive \
     --output=log_o3_profile.txt \
     L1_o3_profile

exit_status="$?"
echo $exit_status > exit_status
exit $exit_status
