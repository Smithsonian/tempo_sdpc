#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
#set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

task_id=$1
run_dir=$2

block_run_subdir=$(printf "block_%03d" $task_id)
cd "${run_dir}/O3PROF/${block_run_subdir}"

export PGSMSG="${SDPC_ROOT}/msgs"
export PGS_PC_INFO_FILE="o3_profile.pcf"

: "${SDPC_O3PROF_TIME_LIMIT:=240}"

srun --ntasks=1 --cpus-per-task=1 \
     --time=$SDPC_O3PROF_TIME_LIMIT \
     --job-name=O3PROF \
     --output=log_o3_profile.txt \
     L1_o3_profile

exit_status="$?"

if test X"$exit_status" != X0 ; then
   printf "\nsrun exit_status=$exit_status\n" >> log_o3_profile.txt
fi

echo $exit_status > exit_status
exit $exit_status
