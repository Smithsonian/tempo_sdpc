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
   echo "${PROGNAME}[$$]: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

log_message()
{
   printf "${PROGNAME}[$$]: $1\n"
}

test -r $rad_path || error_exit "$LINENO: cannot access granule: $rad_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

export PATH="$SDPC_ROOT/bin:$PATH"

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

# Parse the path to the post-INR radiance file
rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

cache_uncompressed_snow_file()
{
   path=$1

   file_sans_gz=$(basename $path .gz)

   cache_dir="$SDPC_RUN_DIR_MASTER/cache/snow"
   if ! test -d $cache_dir ; then
      mkdir -p $cache_dir || error_exit "$LINENO: cannot create directory: $cache_dir"
   fi

   snow_file="$cache_dir/$file_sans_gz"
   if ! test -f "$snow_file" ; then
      gunzip -c $path > "$snow_file" || error_exit "$LINENO: uncompressing file into $cache_dir"
   fi
}

# Generate file list file on master node
irr_file=$(select_irr.py "$rad_path")
snow_file=$(select_ims.py "$rad_path")
cache_uncompressed_snow_file $snow_file
granule_dir=$(dirname "$rad_path")
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
EOF

log_message "start level1b_batch: $SDPC_GRANULE_LABEL"

# Run the post-INR pipeline to prepare for L2 product generation:
sbatch --wait --job-name="L1b" --comment=$SDPC_GRANULE_LABEL \
        --chdir $l1_run_dir \
        --nodes=1-1 --ntasks=8 \
        level1b_batch.sh "${rad_basename}.nc" "$file_list_file"

log_message "finished level1b_batch: $SDPC_GRANULE_LABEL"
