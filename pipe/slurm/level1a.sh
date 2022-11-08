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

# SDPC_NODE_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_NODE_DIR:?SDPC_NODE_DIR not set}"

: "${SDPC_PIPE_ID:?SDPC_PIPE_ID not set}"

# Potential wait times for corresponding Level 0 IRU, SMC, HK
wait_iru_sec="$(config_setting level1a.wait_iru_sec)"
wait_hk_sec="$(config_setting level1a.wait_hk_sec)"

make_iru_only_file_for_inr()
{
   time_interval_file="$1"

   # read time interval from csv file
   OLDIFS=$IFS
   export IFS=','
   read tbeg tend epoch tbeg_utc<"$time_interval_file"
   export IFS=$OLDIFS

   # INR SW will fail if the GOES source directory paths are not found
   tstart="$(date -u --date @$tbeg_utc +%Y-%m-%dT%H:%M:%SZ)"
   prep_inr_goes_source $tstart

   # run L1_inr_prep to generate the INR input file
   export SDPC_NODE_DIR="$SDPC_PIPE_DIR"
   etc_dir="$SDPC_PIPE_DIR/etc"

   echo "run L1_inr_prep: (tbeg,tend)=($tbeg,$tend): $time_interval_file"

   # L1_inr_prep wants file lists in files with standard names.
   # We could edit the config file and change the file list names,
   # or we could run in a subdirectory with a unique name.
   # Creating the subdirectory seems simpler.
   this_dir=$(pwd)
   dir=$(mktemp -d -p $this_dir)
   cd $dir

   select_l0.py --wait $wait_iru_sec --table IRU_L0 --begin $tbeg --end $tend > iru.lis
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

set_dirpath_symlink()
{
   from_dir="$1"
   to_symlink="$2"

   current_path=""
   if test -h "$to_symlink" ; then
      current_path=$(readlink -m "$to_symlink")
   fi

   # create the symlink whether or not $from_dir exists:
   if test x"$current_path" != x"$from_dir" ; then
      ln -nf -s "$from_dir" "$to_symlink"  || error_exit "$LINENO: setting symlink $from_dir -> $to_symlink"
   fi
}

prep_inr_goes_source()
{
   # During normal operations, the GOES imagery source will change only
   # once per day, but the best way to automate this update is to set it here.
   # If the necessary GOES imagery doesn't exist, we proceed with a warning
   # rather than a fatal error, because we don't want to prevent the system
   # from running in an unusual mode (e.g. maybe we don't care about INR
   # in some context, so the lack of imagery is irrelevant).

   # YYYY-MM-DDTHH:mm:ssZ
   tstart="$1"

   yday_subdir="$(TZ='UTC+6' date -d $tstart +%Y/%j)"
   goes_srcdir="${SDPC_ANCILLARY_ROOT}/var/goes/${yday_subdir}"
   if ! test -d "$goes_srcdir" ; then
      echo "WARNING: INR reference GOES imagery not found: $goes_srcdir"
   fi
   target_dir="$SDPC_PIPE_DIR/inr/Staging"

   set_dirpath_symlink $goes_srcdir/east_cmi $target_dir/Right
   set_dirpath_symlink $goes_srcdir/west_cmi $target_dir/Left
}

# Parse the path to the granule file:
granule_basename=$(basename "$granule_path" .nc| sed -e s"/^[.]//")
granule_dir=$(dirname "$granule_path")

iru_file_list=NONE
smc_file_list=NONE
iers_bulletin=NONE

case "${granule_basename}" in
   *_INR_* )
   make_iru_only_file_for_inr $granule_path
   exit 0
   ;;

   *IRR* )
   dark_file_path=$(select_dark.py "$granule_path")
   ntasks=10
   ;;

   *RAD* )
   dark_file_path=$(select_dark.py "$granule_path")
   iers_bulletin=$(select_iers.py "$granule_path")

   iru_file_list="$granule_dir/.${granule_basename}_iru.lis"
   select_l0.py --wait $wait_iru_sec --table IRU_L0 --granule "$granule_path" > $iru_file_list

   smc_file_list="$granule_dir/.${granule_basename}_smc.lis"
   select_l0.py --table SMC_L0 --granule "$granule_path" > $smc_file_list

   tstart=$(global_attribute.py --attr time_coverage_start "$granule_path")
   prep_inr_goes_source $tstart
   ntasks=1
   ;;

   * )
   dark_file_path=NONE
   ntasks=1
   ;;
esac

hk_file_list="$granule_dir/.${granule_basename}_hk.lis"
select_l0.py --wait $wait_hk_sec --table HK_L0 --granule "$granule_path" > $hk_file_list

# Create file-list file
file_list_file="$granule_dir/.${granule_basename}.lis"
cat <<EOF > $file_list_file
granule_path=${granule_path}
dark_file_path=${dark_file_path}
hk_file_list=${hk_file_list}
iru_file_list=${iru_file_list}
smc_file_list=${smc_file_list}
iers_bulletin=${iers_bulletin}
EOF

export SDPC_GRANULE_LABEL="$granule_basename"

# Run the pipeline:
# Time-ordered processing is important:
#  * DRK must finish before the relevant IRR or RAD
#  * RAD time sequence is critical for INR
log_message "submitting sbatch/wait level1a_batch.sh: $SDPC_GRANULE_LABEL"

# Singleton dependency requires a job-name unique to this pipeline.
slurm_logdir="$SDPC_PIPE_DIR/log/level1a/slurm"
jid=$(sbatch --wait --dependency=singleton --parsable \
       --ntasks=$ntasks \
       --job-name="L0:serial:${SDPC_PIPE_ID}" \
       --comment "$SDPC_GRANULE_LABEL" \
       --chdir $run_dir \
       --output "$slurm_logdir/${granule_basename}.level1a_batch-%j.out" \
       level1a_batch.sh "${granule_basename}.nc" "$file_list_file")

log_message "completed: sbatch/wait level1a_batch.sh: $SDPC_GRANULE_LABEL [$jid]"
