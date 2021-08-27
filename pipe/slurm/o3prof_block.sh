#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
#set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

run_dir=$1

block_run_subdir=$(printf "block_%03d" $SLURM_ARRAY_TASK_ID)
cd "${run_dir}/O3PROF/${block_run_subdir}"

export PGSMSG="${SDPC_ROOT}/msgs"
export PGS_PC_INFO_FILE="o3_profile.pcf"

#--exclusive
srun --ntasks=1 --cpus-per-task=1 \
     --output=log_o3_profile.txt \
     L1_o3_profile

exit_status="$?"

if test X"$exit_status" != X0 ; then
   printf "\nsrun exit_status=$exit_status\n" >> log_o3_profile.txt
fi

echo $exit_status > exit_status
exit $exit_status
