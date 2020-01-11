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
   echo "${PROGNAME}: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

test -r $granule_path || error_exit "$LINENO: cannot access granule: $granule_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

export PATH="$SDPC_ROOT/bin:$PATH"

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
   etc_dir="$SDPC_ROOT/etc"

   ephem_file_path=$(filedb -c $SDPC_ROOT/etc/filedb.cfg ephemeris --find --sec $tbeg_utc)

   echo "run L1_inr_prep: (tbeg,tend)=($tbeg,$tend): $time_interval_file"

   # No delay is needed here (any IRU coverage padding extends to earlier times)
   L1_inr_prep -v 1 --Version $SDPC_PROCESSING_VERSION \
       --config ${etc_dir}/l1_inr_prep.cfg \
       --begin $tbeg --end $tend --epoch $epoch \
       --ephemeris ${ephem_file_path}

   echo "L1_inr_prep finished"

   # If L1_inr_prep fails, 'set -e' ensures that the script
   # will exit before we can delete the time interval file.
   # If L1_inr_prep succeeds, its ok to delete the file.
   /bin/rm -f $time_interval_file
}

# Parse the path to the granule file:
granule_basename=$(basename "$granule_path" .nc| sed -e s"/^[.]//")
granule_dir=$(dirname "$granule_path")

case "${granule_basename}" in
   *_INR_* )
   make_iru_only_file_for_inr $granule_path
   exit 0
   ;;

   *IRR* | *RAD* )
   dark_file_path=$(select_dark.py "$granule_path")
   ephem_file_path=$(filedb -c $SDPC_ROOT/etc/filedb.cfg ephemeris --find --header "$granule_path")
   ;;
   * )
   dark_file_path=NONE
   ephem_file_path=NONE
   ;;
esac

# Create file-list file
file_list_file="$granule_dir/.${granule_basename}.lis"
cat <<EOF > $file_list_file
  granule_path=${granule_path}
  dark_file_path=${dark_file_path}
  ephem_file_path=${ephem_file_path}
EOF

export SDPC_GRANULE_LABEL="$granule_basename"

echo "start batch run_L0.sh: $SDPC_GRANULE_LABEL"

# Run the pipeline:
job_l0="L0:$SDPC_GRANULE_LABEL"
sbatch --job-name=$job_l0 \
       --chdir $run_dir \
       run_L0.sh "${granule_basename}.nc" "$file_list_file"
