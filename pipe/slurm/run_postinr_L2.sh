#! /bin/sh -v

set -e
set -u

# Note that the radiance file name be "hidden" (may begin with a ".").

if test $# -ne 3 ; then
  echo "Usage: $0 <rad-path> <run_dir> <product-list>"
  exit 1
fi

rad_path="$1"
run_dir="$2"
product_list="$3"

test -r $rad_path || exit 1
test -d "$SDPC_ROOT" || exit 1

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

export PATH="$SDPC_ROOT/bin:$PATH"

# Parse the path to the post-INR radiance file
rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# construct granule label string for slurm job names
IFS='_' read -r prefix date scan_num granule_num version suffix<<<$rad_basename
scan_num=$(echo $scan_num | sed -e s"/^0*//")
granule_num=$(echo $granule_num | sed -e s"/^0*//")
export SDPC_GRANULE_LABEL="${scan_num}:${granule_num}"

# Run the pipeline:
job_prep_l2="${SDPC_GRANULE_LABEL}:L2-pre"
jid=$(sbatch --parsable --job-name=$job_prep_l2 --chdir $run_dir prep_L2.sh "$rad_path" "${rad_basename}.nc")

job_run_l2="${SDPC_GRANULE_LABEL}:L2"
tarfile_path="$SDPC_RUN_DIR/L2/incoming/${rad_basename}.tar"
sbatch --dependency=afterany:$jid \
       --job-name=$job_run_l2 \
       --chdir $run_dir \
       run_L2.sh "$tarfile_path" "$product_list"
