#! /bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null

# 0. This script is intended to run on a compute node.
#
# 1. Assume this script is started in a writeable directory.
#    The command line arguments are:
#       $1 = granule basename
#       $2 = file containing a list of filenames:
#               granule_path
#               dark_file_path
#
#    The following environment variables are assumed to be set:
#          * SDPC_ROOT, SDPC_RUN_DIR, SDPC_ARCHIVE_DIR
#
# 2. Successful excecution yields the following cases:
#    DRK:  run L0_ccd, archive the result
#    IRR:  run L0_ccd, perform wavelength calibration, archive the result
#    RAD:  run L0_ccd, run L1_inr_prep, archive the (intermediate) result,
#          and put the radiance granule into the INR input cache
#
# 3. On error, a tar file is placed in the repro/L0 directory
#
#---------------------------------------------------------------------

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_RUN_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <target-basename> <list-of-files>"
  exit 1
fi
granule_basename="$1"
file_list_file="$2"

# including this file should define these variables:
#    granule_path
#    dark_file_path
#    hk_file_list
#    iru_file_list
#    smc_file_list
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
etc_dir="$SDPC_RUN_DIR_MASTER/etc"

l0_repro_dir="$SDPC_RUN_DIR_MASTER/repro/L0"
l0_out_dir="$SDPC_RUN_DIR/L0/out"

inr_input_cache="$SDPC_INR_RUN_DIR/Staging/Granules"

# Make a working directory with a local copy of the granule file.
work_dir=$(basename $granule_basename .nc)
/bin/mkdir "$work_dir"
cd $work_dir
/bin/cp "$granule_path" "$granule_basename"
/bin/cp "$hk_file_list" hk.lis
/bin/cp "$file_list_file" "${granule_basename}.lis"
chmod u+w "$granule_basename"

work_dir_tarfile="${work_dir}.tar"

level1_info --dir $granule_basename > archive_subdir

run_dir=$(pwd)
parent_dir=$(dirname "$run_dir")

tar_product_to_dest_dir()
{
   dest_dir=$1
   /bin/mkdir -p $dest_dir

   cd $run_dir
   /bin/rm $granule_basename

   cd $parent_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile $work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile
}

finish()
{
   tar_product_to_dest_dir "$l0_repro_dir"
}
trap finish EXIT ERR

run_l0_ccd()
{
   output_file="$1"
   dark_option="$2"

   /bin/cp ${etc_dir}/l0_ccd.cfg .

   # If no HK files were found, let L0_ccd search the archive for something suitable.
   hk_none=$(grep NONE hk.lis || true)
   if test x"$hk_none" = x"NONE" ; then
      lookup_option=""
   else
      lookup_option="-i @hk.lis"
   fi

   srun --ntasks=1 --output=log_l0_ccd.txt \
   L0_ccd -vv --Version $SDPC_PROCESSING_VERSION $lookup_option \
          --trend trend_params.nc \
          -o $output_file $dark_option \
          $granule_basename
}

. $SDPC_ROOT/bin/wavecal.sh

run_inr_prep()
{
   target_file="$1"

   /bin/cp ${etc_dir}/l1_inr_prep.cfg .

   /bin/cp "$iru_file_list" iru.lis
   /bin/cp "$smc_file_list" smc.lis

   srun --ntasks=1 --output=log_inr_prep.txt \
        L1_inr_prep -v 1 $target_file
}

case "${granule_basename}" in
  *DRK* )
  output_file=$(mkgranule_name -L 1 -p DRK -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file ""
  ;;

  *IRR* )
  # Use irr_type to distinguish IRR (working diffuser) and IRRR (reference diffuser)
  irr_type=$(echo $granule_basename | cut -f2 -d_)
  output_file=$(mkgranule_name -L 1 -p $irr_type -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file "-d $dark_file_path"
  run_wavecal $output_file "$SDPC_IRR_WAVECAL_ARRAY"
  ;;

  *RAD* )
  output_file=$(mkgranule_name -L 1 -p RAD -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file "-d $dark_file_path"
  # run radiance wavelength calibration post-INR, so as not to delay INR
  run_inr_prep $output_file

  rad_tmpfile=$inr_input_cache/.${output_file}
  /bin/cp $output_file $rad_tmpfile
  /bin/mv $rad_tmpfile $inr_input_cache/$output_file
  ;;

  * )
  printf "*** Unsupported granule filename pattern\n"
  exit 1
  ;;
esac

trap - EXIT

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

if test x"$l0_out_dir" != x ; then
   tar_product_to_dest_dir "$l0_out_dir"
   tarfile_path="$l0_out_dir/${work_dir_tarfile}"
   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L1 $tarfile_path
fi

# Assume the initial L0 granule was archived when it was produced,
# so it's ok to delete this copy once the archive.sl process has
# succeeded.
/bin/rm -f "$granule_path" "$file_list_file" "$hk_file_list" "$iru_file_list" "$smc_file_list"
