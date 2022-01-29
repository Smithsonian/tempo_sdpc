#! /bin/sh

# 0. This script is run by cachemon, triggered by the appearance
#    of a newly completed geolocated radiance file in the INR
#    output cache.
#
# 1. The script first runs a batch process (level1b_batch.sh) on a compute
#    node to finish Level 1 processing (post-INR calculations,
#    polarization correction, wavelength calibration), and generate
#    the cloud product.
##
# On error, a tar file is stored in repro/L1
#
#--------------------------------------------------------------------

set -e
set -u

# Note that the radiance file name be "hidden" (may begin with a ".").

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-path> <l1_run_dir>"
  exit 1
fi

rad_path="$1"
l1_run_dir="$2"

PROGNAME="$(basename $0)"
error_exit()
{
   echo "${PROGNAME}[$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

log_message()
{
   printf "${PROGNAME}[$$]: $1\n"
}

test -r $rad_path || error_exit "$LINENO: cannot access granule: $rad_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

# Parse the path to the post-INR radiance file
rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

# Generate file list file on master node
irr_file=$(select_irr.py "$rad_path")
snow_file=$(select_ims.py "$rad_path")
granule_dir=$(dirname "$rad_path")
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
EOF

# Run the post-INR pipeline to prepare for L2 product generation:
slurm_logdir="$SDPC_RUN_DIR_MASTER/log/level1b/slurm"
jid=$(sbatch --job-name="L1b" --parsable \
       --comment=$SDPC_GRANULE_LABEL \
       --chdir $l1_run_dir \
       --output "$slurm_logdir/${rad_basename}.level1b_batch-%j.out" \
       level1b_batch.sh "${rad_basename}.nc" "$file_list_file")

log_message "submitted sbatch $jid: level1b_batch: $SDPC_GRANULE_LABEL"
