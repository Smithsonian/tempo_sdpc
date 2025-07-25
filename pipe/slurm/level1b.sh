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
# Also, "radiance" file may be either TEMPO_RAD_L1 or TEMPO_RADT_L1

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

# When NRT pipeline processing is enabled, create a hard link on the master node
# to trigger NRT processing of this granule
if ! test -f "$SDPC_PIPE_DIR/ctrl/disable-nrt" ; then
   pass2_dir="${SDPC_PIPE_DIR}/stage/granules/inr_output/nrt/inr_pass2"
   if ! test -d $pass2_dir ; then
      mkdir -p $pass2_dir
   fi
   basename_sans_dot=$(basename "$rad_path" | sed -e s"/^[.]//")
   # Don't deliver RADT files to the NRT pipeline - only _RAD_ files.
   case "$basename_sans_dot" in
      *_RAD_* )
        ln $rad_path $pass2_dir/$basename_sans_dot
        ;;
      * )
        ;;
   esac
fi

# SDPC_NODE_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_NODE_DIR:?SDPC_NODE_DIR not set}"

# Parse the path to the post-INR radiance file
rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

# Generate file list file on master node
irr_file=$(select_irr.py --window $SDPC_IRR_SELECT_WINDOW "$rad_path")
snow_file=$(select_ims.py "$rad_path")
granule_dir=$(dirname "$rad_path")

solcal_file_o2o2=""
solcal_file_list="$granule_dir/${rad_basename}.solcal"
# When solar wavelength calilbration caching is enabled,
# generate a list of available cached files.
# Complain if any solcal file is missing when solcal caching is enabled.
if test $SDPC_SOLCAL_CACHE_ENABLE -ne 0 ; then
   # Use cached O2O2 solar wavecal to generate the CLDO4 product.
   solcal_file_o2o2=$(select_irrcal.py --molecule O2O2 $irr_file)
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
         log_message "WARNING: level1b_batch not submitted: $SDPC_GRANULE_LABEL"
         exit 0
      fi
   fi
fi

# Prepare for destriping of trace gas products
destripe_file_list="$granule_dir/${rad_basename}.destripe"
truncate -s 0 $destripe_file_list
if test -n "$SDPC_DESTRIPE_TG" ; then
   _destripe_products=$(echo $SDPC_DESTRIPE_TG | tr , ' ')
   for p in $_destripe_products ; do
       molecule=$(echo $p | tr -d _L2)
       destripe_path=$(select_destripe.py --molecule $molecule ${rad_basename}.nc)
       if test -f "$destripe_path" ; then
          echo $destripe_path >> $destripe_file_list
       fi
   done
fi

# Prepare for CLDO4 destriping
if ! test -f $SDPC_PIPE_DIR/ctrl/disable-destripe-CLDO4 ; then
   destripe_path=$(select_destripe.py --molecule CLDO4 ${rad_basename}.nc)
   if test -f "$destripe_path" ; then
      echo $destripe_path >> $destripe_file_list
   fi
fi

file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
solcal_file_o2o2=${solcal_file_o2o2}
solcal_file_list=${solcal_file_list}
destripe_file_list=${destripe_file_list}
EOF

# If SDPC_RADIANCE_WAVECAL is not set, define it to be ON (non-zero).
# To turn off radiance wavelength calibration, set it to zero.
: "${SDPC_RADIANCE_WAVECAL:=1}"
: "${SDPC_RADIANCE_WAVECAL_NTASKS:=2}"

max_num_tasks=1
if test $SDPC_RADIANCE_WAVECAL -ne 0 ; then
   max_num_tasks=$((2*$SDPC_RADIANCE_WAVECAL_NTASKS))
fi

# Run the post-INR pipeline to prepare for L2 product generation:
slurm_logdir="$SDPC_PIPE_DIR/log/level1b/slurm"
jid=$(sbatch --job-name="L1b" --parsable \
       --comment=$SDPC_GRANULE_LABEL \
       --chdir $l1_run_dir --ntasks=$max_num_tasks \
       --output "$slurm_logdir/${rad_basename}.level1b_batch-%j.out" \
       level1b_batch.sh "${rad_basename}.nc" "$file_list_file")

log_message "submitted sbatch $jid: level1b_batch: $SDPC_GRANULE_LABEL"
