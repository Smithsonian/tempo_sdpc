#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
# set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If SDPC_O3PROF_MODE is not set, define it
: "${SDPC_O3PROF_MODE:=UVVIS}"

# If USE_SYNTHETIC_MET_DATA is not set, define it to be OFF
# To use synthetic met data, set it to the relevant file path
: "${USE_SYNTHETIC_MET_DATA:=OFF}"

# 1. Processing will run in the subdirectory provided on the command line,
#    which already contains all necessary inputs.
# 2. Processing ultimately stores all results in a tar file in an
#    appropriate destination directory.
# 3. When processing ends, remove the processing directory.

mode="$1"
host_spec="$2"

case $mode in
  init | cleanup)
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
l2_repro_dir="$SDPC_RUN_DIR_MASTER/repro/L2"
etc_dir="$SDPC_RUN_DIR_MASTER/etc"

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
lev2_file_fmt=$(mkgranule_name -L 2 -p %s -v $SDPC_PROCESSING_VERSION ${rad_basename}.nc)
lev2_base_fmt=$(basename $lev2_file_fmt .nc)

out_basename=$(printf "$lev2_base_fmt" O3PROF)

tar_product_to_dest_dir()
{
   dest_dir=$1

   cd $run_dir

   if ! test -d $work_dir ; then
      return
   fi

   /bin/rm $work_dir/${rad_basename}.nc
   /bin/rm $work_dir/${irr_basename}.nc
   /bin/rm $work_dir/$cld_file

   granule_dir=$(basename $run_dir)

   # We want every product tar file to contain a copy of the archive_subdir file,
   # and we want tar to delete all the files it collects.
   # To avoid deleting the archive_subdir file we create a temporary directory,
   # copy the archive_subdir file over there, and then create the tar file.

   this_host=$(uname -n | cut -d. -f1)

   cd $parent_dir
   tmp_dir="$(mktemp -d tmpo3p.XXXXXX)"
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
   dirs="$(find . -maxdepth 1 -type d -name 'block_*')"
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
        /bin/rm -f pge_input_basenames.lis ${rad_basename}.lis
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

define_met_files()
{
met_file1="notused"
met_dir1=""
met_file2="notused"
met_dir2=""
if ! test x"$USE_SYNTHETIC_MET_DATA" = x"OFF" ; then
      met_file1=$(basename $USE_SYNTHETIC_MET_DATA)
      met_dir1=$(dirname $USE_SYNTHETIC_MET_DATA)
fi
}

config_subdir()
{
   subdir_name=$1

   case "$SDPC_O3PROF_MODE" in
     UV)
       control_file="${etc_dir}/o3_profile/default_main_control_uv.inp"
       profoz_file="${etc_dir}/o3_profile/profoz_uv.inp"
       ;;

     UVVIS)
       control_file="${etc_dir}/o3_profile/default_main_control_uvvis.inp"
       profoz_file="${etc_dir}/o3_profile/profoz_uvvis.inp"
       ;;

     *)
     echo "*** $0: unsupported mode SDPC_O3PROF_MODE = $SDPC_O3PROF_MODE"
     exit 1
     ;;
   esac

   # Copy control files to product directory
   /bin/cp $control_file $subdir_name
   /bin/cp $profoz_file $subdir_name

   etc_dir="$SDPC_RUN_DIR_MASTER/etc"

   product_dir=.
   spectra_dir=.
   cloud_dir=.
   refdata_dir="$SDPC_RUN_DIR/refdata/o3_profile"

   pcf_file="$product_dir/o3_profile.pcf"

   # file names
   radiance_file="../${rad_basename}.nc"
   irradiance_file="../${irr_basename}.nc"
   cloud_file="../${cld_file}"
   product_file="${out_basename}.nc"
   control_file_basename=$(basename $control_file)
   profoz_file_basename=$(basename $profoz_file)

   define_met_files

   # Read the block parameters:
   read xt_beg xt_end <"$subdir_name/block.txt"
   # processing code uses 1-based numbering, so add 1:
   xt_beg=$(($xt_beg+1))
   xt_end=$(($xt_end+1))

   template_pcf="$etc_dir/o3_profile/o3_profile.pcf.in"
# Edit the PCF file template:
   sed \
    -e s,@refdata_dir@,$refdata_dir,g \
    -e s,@spectra_dir@,$spectra_dir,g \
    -e s,@cloud_dir@,$cloud_dir,g \
    -e s,@product_dir@,$product_dir,g \
    -e s,@etc_dir@,$etc_dir,g \
    -e s,@radiance_file@,$radiance_file,g \
    -e s,@irradiance_file@,$irradiance_file,g \
    -e s,@cloud_file@,$cloud_file,g \
    -e s,@product_file@,$product_file,g \
    -e s,@control_file@,$control_file_basename,g \
    -e s,@profoz_file@,$profoz_file_basename,g \
    -e s,@met_dir1@,$met_dir1,g \
    -e s,@met_file1@,$met_file1,g \
    -e s,@met_dir2@,$met_dir2,g \
    -e s,@met_file2@,$met_file2,g \
    -e s,@line_sample_extent@,-1\ -1\ $xt_beg\ $xt_end,g \
    -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
    $template_pcf > $subdir_name/$pcf_file
}

# The string $host_spec has the form k-N, where 0 <= k < N
# indicating that the calculation is being spread over N hosts,
# and this particular host is the kth one.

array_bounds=$(o3p_partition.sl -m . "$host_spec")

block_dirs="$(find . -maxdepth 1 -type d -name 'block_*')"

for subdir in $block_dirs; do
   config_subdir $subdir
done

trap - EXIT
# Return the array bounds:
# (this is the intended output -- don't delete this!)
echo "$array_bounds"
