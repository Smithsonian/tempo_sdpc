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
export SDPC_GRANULE_LABEL="${rad_basename}"

# FIXME - in operations, a cron job should handle this DB update, e.g. 1x per day
# For testing, we'll update it here, because the irr file may be newly created
filedb -c $SDPC_ROOT/etc/filedb.cfg tempo:irr --update

# Generate file list file on master node
irr_file=$(filedb -c $SDPC_ROOT/etc/filedb.cfg tempo:irr --find --header "$rad_path")
snow_file=$(filedb -c $SDPC_ROOT/etc/filedb.cfg snow --find --header "$rad_path")
met_file_path=$(filedb -c $SDPC_ROOT/etc/filedb.cfg met --find --header "$rad_path")
granule_dir=$(dirname "$rad_path")
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
met_file_path=${met_file_path}
EOF

# Run the pipeline:
job_prep_l2="L2-pre:${SDPC_GRANULE_LABEL}"
jid=$(sbatch --parsable --job-name=$job_prep_l2 --chdir $run_dir prep_L2.sh "${rad_basename}.nc" "$file_list_file")

job_run_l2="L2:${SDPC_GRANULE_LABEL}"
tarfile_path="$SDPC_RUN_DIR/L2/incoming/${rad_basename}.tar"
sbatch --dependency=afterany:$jid \
       --job-name=$job_run_l2 \
       --chdir $run_dir \
       run_L2.sh "$tarfile_path" "$product_list"
