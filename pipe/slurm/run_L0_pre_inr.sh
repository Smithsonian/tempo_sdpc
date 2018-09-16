#! /bin/sh -v

set -e
set -u

# Note that the input file name be "hidden" (may begin with a ".").

if test $# -ne 2 ; then
  echo "Usage: $0 <granule-path> <run_dir>"
  exit 1
fi

granule_path="$1"
run_dir="$2"

test -r $granule_path || exit 1
test -d "$SDPC_ROOT" || exit 1

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

export PATH="$SDPC_ROOT/bin:$PATH"

# Parse the path to the granule file
granule_basename=$(basename "$granule_path" .nc| sed -e s"/^[.]//")

export SDPC_GRANULE_LABEL="$granule_basename"

# Run the pipeline:
job_l0="L0:$SDPC_GRANULE_LABEL"
sbatch --job-name=$job_l0 \
       --chdir $run_dir \
       run_L0.sh "$granule_path" "${granule_basename}.nc"
