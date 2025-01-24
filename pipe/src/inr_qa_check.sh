#! /bin/sh

# Perform INR quality check:
#  - Create a temporary working directory in $output_dir
#  - Populate the working directory with subdirectories containing
#    symlinks to radiance data and GOES imagery for a single day
#  - Run OTS code to generate diagnostic plots
#  - Create a tar file containing the diagnostic output
#  - Delete the working directory

set -u

PROGNAME="$(basename $0)"

if test $# -ne 2 ; then
   echo "Usage:  $PROGNAME YYYY-MM-DD <output-dir>"
   exit 0
fi

check_cmd_exists()
{
   cmd="$1"
   if ! command -v $cmd &> /dev/null ; then
      echo "*** Error: $cmd not found"
      exit 1
   fi
}

check_cmd_exists tempo_inr_quality.sh
check_cmd_exists inr_qa_setup.py

# Date must have the form: YYYY-MM-DD
date_ymd="$1"
output_dir="$2"

: "${SDPC_ARCHIVE_DBFILE:?SDPC_ARCHIVE_DBFILE not set}"
: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

error_exit()
{
   echo "${PROGNAME}: ERROR: ${1:-'Unknown Error'}" 1>&2
   exit 1
}
trap error_exit ERR

# If necessary, create the output directory
if ! test -d "$output_dir" ; then
   mkdir -p "$output_dir"
fi

# Create temporary working directory
work_dir="$(mktemp -p $output_dir -d)"
mkdir -p "$work_dir"

# Populate the working directory with symlinks to input data
inr_qa_setup.py --dbfile "$SDPC_ARCHIVE_DBFILE" --dir "$work_dir" $date_ymd || error_exit "inr_qa_setup.py failed"

# Unset DISPLAY to eliminate irrelevant complaints
# from OTS SW about about port 6011 connections refused
unset DISPLAY

# Generate diagnostic plots
tempo_inr_quality.sh "$work_dir/config.txt" || error_exit "tempo_inr_quality.sh failed"

# Name the tar file using the satellite day number extracted from the first radiance file
first_radiance_path=$(find $work_dir/radiances -maxdepth 1 -type l | sort | head -n 1)
sat_day=$(level1_info --localday $first_radiance_path | tr -d D)
tar_basename="tempo_inrq_d${sat_day}"
tar_path="${output_dir}/${tar_basename}.tar"

# Collect only the diagnostic output in a tar file (does not compress well so dont bother)
tar cf "$tar_path" -C "$work_dir" output --transform="s,output,${tar_basename}," || error_exit "tar failed"

# Delete entire working directory after successful completion
if test -d "$work_dir" ; then
   /bin/rm -rf "$work_dir"
fi
