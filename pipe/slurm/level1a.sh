#! /bin/sh

# 0. This script is run by cachemon, triggered by the arrival of
#    each new granule of exposure records, whether dark, irradiance,
#    or radiance.  This script may also be triggered by the arrival
#    of a special "INR" file containing a time interval for which IRU
#    coverage should be forwarded to the INR subsystem.
#
# 1. When triggered on behalf of a granule of Level 0 exposure records,
#    this script constructs a list of data files needed to support
#    processing and then submits a slurm batch job to process the
#    Level 0 data (in a subdirectory of the directory specified on the
#    command line).
#    When triggered on behalf of an "INR" file, this script runs
#    L1_inr_prep to collect the relevant IRU, SMC, and ephemeris time
#    series in a special "telemetry only" granule for input to INR.

set -e
set -u

# Note that the input file name be "hidden" (may begin with a ".").

if test $# -ne 2 ; then
  echo "Usage: $0 <granule-path> <run_dir>"
  exit 1
fi

granule_path="$1"
run_dir="$2"

PROGNAME="$(basename $0)"
error_exit()
{
   echo "${PROGNAME}: ERROR: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

log_message()
{
   printf "${PROGNAME}[$$]: $1\n"
}

test -r $granule_path || error_exit "$LINENO: cannot access granule: $granule_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

make_iru_only_file_for_inr()
{
   time_interval_file="$1"

   # read time interval from csv file
   OLDIFS=$IFS
   export IFS=','
   read tbeg tend epoch tbeg_utc<"$time_interval_file"
   export IFS=$OLDIFS

   # run L1_inr_prep to generate the INR input file
   export SDPC_RUN_DIR="$SDPC_RUN_DIR_MASTER"
   etc_dir="$SDPC_RUN_DIR_MASTER/etc"

   echo "run L1_inr_prep: (tbeg,tend)=($tbeg,$tend): $time_interval_file"

   # L1_inr_prep wants file lists in files with standard names.
   # We could edit the config file and change the file list names,
   # or we could run in a subdirectory with a unique name.
   # Creating the subdirectory seems simpler.
   this_dir=$(pwd)
   dir=$(mktemp -d -p $this_dir)
   cd $dir

   select_l0.py --wait 120 --table IRU_L0 --begin $tbeg --end $tend > iru.lis
   select_l0.py --table SMC_L0 --begin $tbeg --end $tend > smc.lis
   select_l0.py --table HK_L0  --begin $tbeg --end $tend > hk.lis

   # Output goes to inr input cache.
   L1_inr_prep -v 1 --Version $SDPC_PROCESSING_VERSION \
       --config ${etc_dir}/l1_inr_prep.cfg \
       --begin $tbeg --end $tend --epoch $epoch

   echo "L1_inr_prep finished"

   # If L1_inr_prep fails, 'set -e' ensures that the script
   # will exit before we can delete the time interval file.
   # If L1_inr_prep succeeds, its ok to delete the file.
   /bin/rm -f $time_interval_file iru.lis smc.lis hk.lis
   cd $this_dir
   /bin/rmdir $dir
}

prep_inr_goes_source()
{
   # During normal operations, the GOES imagery source will change only
   # once per day, but the best way to automate this update is to set it
   # here, using the radiance file measurement date. While we could check
   # the existing symlinks before changing them, the code is simplest
   # if we set the symlinks every time.  Re-setting the links for each new
   # radiance file also ensures a quick recovery in case something happens
   # to the links during the course of an operational day.
   # If the necessary GOES imagery doesn't exist, we proceed with a warning
   # rather than a fatal error, because we don't want to prevent the system
   # from running in an unusual mode (e.g. maybe we don't care about INR
   # in some context, so the lack of imagery is irrelevant).

   tstart=$(radiance_attribute.py --attr time_coverage_start "$granule_path")
   yday_subdir="$(TZ='UTC+6' date -d $tstart +%Y/%j)"
   goes_srcdir="${SDPC_ANCILLARY_ROOT}/goes/${yday_subdir}"

   if test -d "$goes_srcdir" ; then
      target_dir="$SDPC_RUN_DIR_INR/Staging"
      ln -nf -s $goes_srcdir/g16_cmi $target_dir/Right || error_exit "$LINENO: setting GOES-East source"
      ln -nf -s $goes_srcdir/g17_cmi $target_dir/Left || error_exit "$LINENO: setting GOES-West source"
   else
      echo "WARNING: INR reference GOES imagery not found: $goes_srcdir"
   fi
}

# Parse the path to the granule file:
granule_basename=$(basename "$granule_path" .nc| sed -e s"/^[.]//")
granule_dir=$(dirname "$granule_path")

iru_file_list=NONE
smc_file_list=NONE

case "${granule_basename}" in
   *_INR_* )
   make_iru_only_file_for_inr $granule_path
   exit 0
   ;;

   *IRR* )
   dark_file_path=$(select_dark.py "$granule_path")
   ;;

   *RAD* )
   dark_file_path=$(select_dark.py "$granule_path")

   iru_file_list="$granule_dir/.${granule_basename}_iru.lis"
   select_l0.py --table IRU_L0 --granule "$granule_path" > $iru_file_list

   smc_file_list="$granule_dir/.${granule_basename}_smc.lis"
   select_l0.py --table SMC_L0 --granule "$granule_path" > $smc_file_list

   prep_inr_goes_source
   ;;

   * )
   dark_file_path=NONE
   ;;
esac

hk_file_list="$granule_dir/.${granule_basename}_hk.lis"
select_l0.py --table HK_L0 --granule "$granule_path" > $hk_file_list

# Create file-list file
file_list_file="$granule_dir/.${granule_basename}.lis"
cat <<EOF > $file_list_file
granule_path=${granule_path}
dark_file_path=${dark_file_path}
hk_file_list=${hk_file_list}
iru_file_list=${iru_file_list}
smc_file_list=${smc_file_list}
EOF

export SDPC_GRANULE_LABEL="$granule_basename"

# Run the pipeline:
# Time-ordered processing is important:
#  * DRK must finish before the relevant IRR or RAD
#  * RAD time sequence is critical for INR
log_message "submitting sbatch/wait level1a_batch.sh: $SDPC_GRANULE_LABEL"

slurm_logdir="$SDPC_RUN_DIR_MASTER/log/level1a/slurm"
jid=$(sbatch --wait --dependency=singleton --parsable \
       --job-name="L0:serial" \
       --comment "$SDPC_GRANULE_LABEL" \
       --chdir $run_dir \
       --output "$slurm_logdir/${granule_basename}.level1a_batch-%j.out" \
       level1a_batch.sh "${granule_basename}.nc" "$file_list_file")

log_message "completed: sbatch/wait level1a_batch.sh: $SDPC_GRANULE_LABEL [$jid]"
