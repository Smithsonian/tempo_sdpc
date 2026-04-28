#! /bin/sh

# 0. This script is run by cachemon, triggered by the appearance
#    of a newly completed first-pass geolocated radiance file in
#    the INR output cache.
#
# 1. The script first runs a batch process (level1b_nrt1_batch.sh) on a compute
#    node to perform the first step of CLDO4 retrieval
##
# On error, a tar file is stored in repro/L1
#
#--------------------------------------------------------------------

set -u

# Note that the radiance file name be "hidden" (may begin with a ".").

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-path> <l1_run_dir>"
  exit 1
fi

orig_rad_path="$1"
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

test -r $orig_rad_path || error_exit "$LINENO: cannot access granule: $orig_rad_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# SDPC_NODE_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_NODE_DIR:?SDPC_NODE_DIR not set}"

# Parse the path to the geolocated radiance file
orig_rad_basename=$(basename "$orig_rad_path" .nc| sed -e s"/.NavigatedResult$//" -e s"/^[.]//")

# create NRT file basename
rad_basename=$(echo $orig_rad_basename | sed -e s"/_V/_NRT_V/")
# Because the version numbers of the baseline and NRT products may differ,
# the version number in $rad_basename may be incorrect. Fix the filename here
# (the version_id attribute in the file gets fixed elsewhere).
current_version_id=$(global_attribute.py --attr version_id $orig_rad_path)
if test $current_version_id -ne $SDPC_NRT_PROCESSING_VERSION ; then
   current_version_label=$(printf "V%02d" $current_version_id)
   new_version_label=$(printf "V%02d" $SDPC_NRT_PROCESSING_VERSION)
   rad_basename=$(echo $rad_basename | sed -e "s/${current_version_label}/${new_version_label}/")
fi

# Create a hard link to preserve a copy of the input radiance file
stage_dir="$SDPC_PIPE_DIR/stage/granules/inr_output/nrt/inr_pass1"
if ! test -d $stage_dir ; then
   mkdir -p $stage_dir
fi
rad_path="$stage_dir/${rad_basename}.nc"
ln $orig_rad_path $rad_path

# Create a directory to receive the result tar file
cldo4_input_dir="$SDPC_PIPE_DIR/stage/granules/cldo4_input_nrt"
if ! test -d $cldo4_input_dir ; then
   mkdir -p $cldo4_input_dir
fi

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

# Prepare file list file on master node
irr_file=$(select_irr.py --window $SDPC_IRR_SELECT_WINDOW "$rad_path")
snow_file=$(select_ims.py "$rad_path")
granule_dir=$(dirname "$rad_path")

# When solar wavelength calilbration caching is enabled,
# we'll use the O2O2 file for generating the CLDO4 product.
solcal_file_o2o2=""
if test $SDPC_SOLCAL_CACHE_ENABLE -ne 0 ; then
   solcal_file_o2o2=$(select_irrcal.py --molecule O2O2 $irr_file)
fi

# Write file list file on master node
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
orig_rad_path=${orig_rad_path}
solcal_file_o2o2=${solcal_file_o2o2}
EOF

# Turn off radiance wavelength calibration:
: "${SDPC_RADIANCE_WAVECAL:=0}"

# Use the first-pass radiance file to complete the first step of CLDO4 retrieval
slurm_logdir="$SDPC_PIPE_DIR/log/level1b_nrt1/slurm"
jid=$(sbatch --job-name="nL1b.1" --parsable --partition="$SDPC_NRT_PARTITION" \
       --comment=$SDPC_GRANULE_LABEL --quiet \
       --chdir $l1_run_dir --ntasks=1 \
       --output "$slurm_logdir/${rad_basename}.level1b_nrt1_batch-%j.out" \
       level1b_nrt1_batch.sh "${rad_basename}.nc" "$file_list_file")

log_message "submitted sbatch $jid: level1b_nrt1_batch: $SDPC_GRANULE_LABEL"
