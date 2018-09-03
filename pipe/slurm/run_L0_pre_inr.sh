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

dark_handler(){
  export SDPC_GRANULE_LABEL="drk0"
}

irradiance_handler(){
  export SDPC_GRANULE_LABEL="irr0"
}

radiance_handler(){
  # construct granule label string for slurm job names
  IFS='_' read -r prefix date scan_num granule_num version suffix<<<$granule_basename
  scan_num=$(echo $scan_num | sed -e s"/^0*//")
  granule_num=$(echo $granule_num | sed -e s"/^0*//")
  export SDPC_GRANULE_LABEL="rad0:${scan_num}:${granule_num}"
}

case "$granule_basename" in
  *drk0* )
  dark_handler
  ;;

  *irr0* )
  irradiance_handler
  ;;

  *rad0* )
  radiance_handler
  ;;

  * )
  echo "*** Error: unsupported filename pattern: $granule_basename"
  exit 1
  ;;
esac

# Run the pipeline:
job_l0="$SDPC_GRANULE_LABEL"
sbatch --job-name=$job_l0 \
       --chdir $run_dir \
       run_L0.sh "$granule_path" "${granule_basename}.nc"
