#! /bin/sh
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1

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
#          * SDPC_ROOT, SDPC_NODE_DIR, SDPC_ARCHIVE_DIR
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
test -d $SDPC_NODE_DIR || exit 1
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
#    irr_file_path
#    hk_file_list
#    iru_file_list
#    smc_file_list
#    iers_bulletin
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
etc_dir="$SDPC_PIPE_DIR/etc"

l0_repro_dir="$SDPC_PIPE_DIR/repro/L0"
l0_out_dir="$SDPC_NODE_DIR/L0/out"

inr_input_cache="$SDPC_PIPE_DIR/inr/Staging/Granules"

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

   case "${granule_basename}" in
      *DRK* | *RADT* )
          config_file="l0_ccd_drk.cfg"
          ;;
      * )
          config_file="l0_ccd.cfg"
          ;;
   esac
   /bin/cp ${etc_dir}/${config_file} .

   if test -n "$irr_file_path" ; then
      solar_option="--solar $irr_file_path"
   else
      solar_option=""
   fi

   # If no HK files were found, let L0_ccd search the archive for something suitable.
   hk_none=$(grep NONE hk.lis || true)
   if test x"$hk_none" = x"NONE" ; then
      lookup_option=""
   else
      lookup_option="-i @hk.lis"
   fi

   srun --ntasks=1 --output=log_l0_ccd.txt \
   L0_ccd -vv --Version $SDPC_PROCESSING_VERSION \
          --config $config_file $lookup_option $solar_option \
          --trend trend_params.nc $dark_option \
          -o $output_file \
          $granule_basename
}

run_inr_prep()
{
   target_file="$1"

   /bin/cp ${etc_dir}/l1_inr_prep.cfg .

   /bin/cp "$iru_file_list" iru.lis
   /bin/cp "$smc_file_list" smc.lis

   if test -f "$iers_bulletin" ; then
      iers_opt="-i $iers_bulletin"
   else
      iers_opt=""
   fi

   srun --ntasks=1 --output=log_inr_prep.txt \
        L1_inr_prep -v 1 $iers_opt $target_file
}

cache_tracegas_solcal()
{
   irr_file=$1

   # Generate these files only when the working diffuser was used
   case "$irr_file" in
      *IRRR* )
          return
          ;;
      * )
          ;;
   esac

   product_list="$(echo $SDPC_SOLCAL_CACHE_PRODUCTS | tr -s , ' ')"
   if test x"$product_list" = x ; then
      return
   fi

   irr_basename=$(basename $irr_file .nc)

   # Run a background process for each trace gas product:
   slurm_logdir="$SDPC_PIPE_DIR/log/level1a/slurm"
   for molecule in $product_list ; do
     tracegas_log="$slurm_logdir/${irr_basename}.tracegas-$molecule-${SLURM_JOB_ID}.out"
     tracegas_solcal_cache.sh $molecule $irr_file > $tracegas_log 2>&1 &
   done

   # wait for all background jobs to exit
   wait
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
  wavecal.sh $output_file 5
  if test $SDPC_SOLCAL_CACHE_ENABLE -ne 0 ; then
     cache_tracegas_solcal $output_file
  fi
  ;;

  *RAD* )
  # Use rad_type to distinguish RAD (radiance) and RADT (twilight scan)
  rad_type=$(echo $granule_basename | cut -f2 -d_)
  output_file=$(mkgranule_name -L 1 -p $rad_type -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file "-d $dark_file_path"
  # run radiance wavelength calibration post-INR, so as not to delay INR
  run_inr_prep $output_file

  # Both RAD and RADT are delivered to the INR input cache
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
