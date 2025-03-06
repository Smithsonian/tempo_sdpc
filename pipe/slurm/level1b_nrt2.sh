#! /bin/sh

# 0. This script is run by cachemon, triggered by the appearance
#    of a newly completed first-pass geolocated radiance file in
#    the INR output cache.
#
# 1. The script first runs a batch process (level1b_nrt2_batch.sh) on a compute
#    node to perform the last step of CLDO4 retrieval, and create a tar file
#    to be used as input for NRT Level 2 processing
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

# SDPC_NODE_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_NODE_DIR:?SDPC_NODE_DIR not set}"

# Parse the path to the geolocated radiance file
orig_rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# create NRT file basename
rad_basename=$(echo $orig_rad_basename | sed -e s"/_V/_NRT_V/")
# Because the version numbers of the baseline and NRT products may differ,
# the version number in $rad_basename may be incorrect. Fix the filename here
# (the version_id attribute in the file gets fixed elsewhere).
current_version_id=$(global_attribute.py --attr version_id $rad_path)
if test $current_version_id -ne $SDPC_NRT_PROCESSING_VERSION ; then
   current_version_label=$(printf "V%02d" $current_version_id)
   new_version_label=$(printf "V%02d" $SDPC_NRT_PROCESSING_VERSION)
   rad_basename=$(echo $rad_basename | sed -e "s/${current_version_label}/${new_version_label}/")
fi

# If necessary, create a directory to receive the result tar notice file
l2_incoming_nrt="$SDPC_PIPE_DIR/stage/granules/level2_input_nrt"
if ! test -d $l2_incoming_nrt ; then
   mkdir -p $l2_incoming_nrt
fi

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

irr_select_window=$(config_setting level1b.irr_select_window)

# Prepare file list file on master node
irr_file=$(select_irr.py --window $irr_select_window "$rad_path")
snow_file=$(select_ims.py "$rad_path")
granule_dir=$(dirname "$rad_path")

solcal_file_list="$granule_dir/${rad_basename}.solcal"
# When solar wavelength calilbration caching is enabled,
# generate a list of available cached files.
# Complain if any solcal file is missing when solcal caching is enabled.
if test $SDPC_SOLCAL_CACHE_ENABLE -ne 0 ; then
   product_list="$(echo $SDPC_SOLCAL_CACHE_PRODUCTS | tr -s , ' ')"
   if test x"$product_list" != x ; then
      truncate -s 0 $solcal_file_list
      missing_solcal_file=0
      for molecule in $product_list ; do
          solcal_path=$(select_irrcal.py --molecule $molecule $irr_file)
          if test -f "$solcal_path" ; then
             echo $solcal_path >> $solcal_file_list
          else
             log_message "WARNING: no $molecule solcal file for: $irr_file"
             missing_solcal_file=1
          fi
      done
      if test $missing_solcal_file -ne 0 ; then
         log_message "WARNING: level1b_nrt2_batch not submitted: $SDPC_GRANULE_LABEL"
         exit 0
      fi
   fi
fi

# Write file list file on master node
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
solcal_file_list=${solcal_file_list}
EOF

# Turn off radiance wavelength calibration:
: "${SDPC_RADIANCE_WAVECAL:=0}"

# Use the second-pass radiance file to finish prep for L2 retrievals
slurm_logdir="$SDPC_PIPE_DIR/log/level1b_nrt2/slurm"
jid=$(sbatch --job-name="nL1b.2" --parsable --partition="$SDPC_NRT_PARTITION" \
       --comment=$SDPC_GRANULE_LABEL \
       --chdir $l1_run_dir --ntasks=1 \
       --output "$slurm_logdir/${rad_basename}.level1b_nrt2_batch-%j.out" \
       level1b_nrt2_batch.sh "${rad_basename}.nc" "$file_list_file")

log_message "submitted sbatch $jid: level1b_nrt2_batch: $SDPC_GRANULE_LABEL"
