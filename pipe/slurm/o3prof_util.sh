#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# 1. Processing will run in the subdirectory provided on the command line,
#    which already contains all necessary inputs.
# 2. Processing ultimately stores all results in a tar file in an
#    appropriate destination directory.
# 3. When processing ends, remove the processing directory.

mode="$1"
host_spec="$2"
block_range_file=""

case $mode in
  init)
    block_range_file="$3"
    ;;
  cleanup)
    ;;
  *)
    echo "*** $0: unsupported mode = $mode"
    exit 1
    ;;
esac

work_dir="O3PROF"
run_dir=$(pwd)
parent_dir=$(dirname $run_dir)
cd $work_dir

l2_out_dir="$SDPC_RUN_DIR/L2/out"
l2_repro_dir="$SDPC_RUN_DIR/L2/repro"
etc_dir="$SDPC_ROOT/etc"

# get input file names
. ./pge_input_basenames.lis
rad_file=$RAD
irr_file=$IRR
cld_file=$CLD

rad_basename=$(basename $rad_file .nc)
irr_basename=$(basename $irr_file .nc)

work_dir_tarfile="${rad_basename}.${work_dir}_${host_spec}.tar"

# Define product file name template
#
lev2_file_fmt=$(mkgranule_name -L 2 -p %s ${rad_basename}.nc)
lev2_base_fmt=$(basename $lev2_file_fmt .nc)

out_basename=$(printf "$lev2_base_fmt" O3PROF)

tar_product_to_dest_dir()
{
   dest_dir=$1

   cd $run_dir

   /bin/rm $work_dir/${rad_basename}.nc
   /bin/rm $work_dir/${irr_basename}.nc
   /bin/rm $work_dir/$cld_file

   granule_dir=$(basename $run_dir)

   # We want every product tar file to contain a copy of the archive_subdir file,
   # and we want tar to delete all the files it collects.
   # To avoid deleting the archive_subdir file we create a temporary directory,
   # copy the archive_subdir file over there, and then create the tar file.

   cd $parent_dir
   tmp_dir="$(mktemp -d tmp.XXXXXX)"
   tmp_granule_dir="$tmp_dir/$granule_dir"
   /bin/mkdir -p $tmp_granule_dir
   /bin/cp $granule_dir/archive_subdir $tmp_granule_dir
   /bin/mv $granule_dir/$work_dir $tmp_granule_dir

   cd $tmp_dir

   /bin/mkdir -p $dest_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile \
         $granule_dir/archive_subdir $granule_dir/$work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile

   cd $parent_dir
   /bin/rmdir $tmp_granule_dir $tmp_dir
}

decide_cleanup_dest_dir()
{
   dirs="$(ls -d block_*)"
   cleanup_dest_dir="$l2_out_dir"
   for d in $dirs ; do
     status_file="$d/exit_status"
     if ! test -r $status_file ; then
        cleanup_dest_dir="$l2_repro_dir"
        break
     fi
     s=$(cat $status_file)
     if test x"$s" != x0 ; then
        cleanup_dest_dir="$l2_repro_dir"
        break
     fi
   done
}

case $mode in
  cleanup)
    decide_cleanup_dest_dir

    if test "$cleanup_dest_dir" != "$l2_out_dir" ; then
        tar_product_to_dest_dir "$cleanup_dest_dir"
    else
        /bin/rm -f pge_input_basenames.lis blocks ${rad_basename}.lis
        tar_product_to_dest_dir "$cleanup_dest_dir"
        tarfile_path="$l2_out_dir/$work_dir_tarfile"
        archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
    fi
    exit "$?"
    ;;
  *)
    ;;
esac

#--------------------------------------------------
#  The code below is used only when mode != cleanup
#--------------------------------------------------

finish()
{
   tar_product_to_dest_dir "$l2_repro_dir"
}
trap finish EXIT ERR

error_exit(){
  echo $1 >&2
  exit 1
}

config_subdir()
{
   subdir_name=$1

   # Copy control files to product directory
   control_file="${etc_dir}/o3_profile/default_main_control.inp"
   /bin/cp $control_file $subdir_name
   profoz_file="${etc_dir}/o3_profile/default_profoz.inp"
   /bin/cp $profoz_file $subdir_name/profoz.inp

   product_dir=.
   spectra_dir=.
   cloud_dir=.
   refdata_dir="$SDPC_RUN_DIR/refdata/o3_profile_jbak"

   pcf_file="$product_dir/o3_profile.pcf"

   # file names
   radiance_file="../${rad_basename}.nc"
   irradiance_file="../${irr_basename}.nc"
   cloud_file="../${cld_file}"
   product_file="${out_basename}.nc"
   control_file_basename=$(basename $control_file)

   # Read the block parameters:
   read xt_beg xt_end <"$subdir_name/block.txt"
   # processing code uses 1-based numbering, so add 1:
   xt_beg=$(($xt_beg+1))
   xt_end=$(($xt_end+1))

   template_pcf="$etc_dir/o3_profile/default.pcf.in"
# Edit the PCF file template:
   sed \
    -e s,@refdata_dir@,$refdata_dir,g \
    -e s,@spectra_dir@,$spectra_dir,g \
    -e s,@cloud_dir@,$cloud_dir,g \
    -e s,@product_dir@,$product_dir,g \
    -e s,@radiance_file@,$radiance_file,g \
    -e s,@irradiance_file@,$irradiance_file,g \
    -e s,@cloud_file@,$cloud_file,g \
    -e s,@product_file@,$product_file,g \
    -e s,@control_file@,$control_file_basename,g \
    -e s,@line_sample_extent@,-1\ -1\ $xt_beg\ $xt_end,g \
    -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
    $template_pcf > $subdir_name/$pcf_file
}

# The string $host_spec has the form k-N, where 0 <= k < N
# indicating that the calculation is being spread over N hosts,
# and this particular host is the kth one.

o3p_partition.sl -m . -b $block_range_file "$host_spec"

block_dirs=$(ls -d block_*)

for subdir in $block_dirs; do
   config_subdir $subdir
done

trap - EXIT
