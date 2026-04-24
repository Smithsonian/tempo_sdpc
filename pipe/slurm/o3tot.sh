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

work_dir="O3TOT"

l2_out_dir="$SDPC_NODE_DIR/L2/out"
l2_repro_dir="$SDPC_PIPE_DIR/repro/L2"

run_dir=$(pwd)
parent_dir=$(dirname $run_dir)
cd $work_dir

# get input file names
. ./pge_input_basenames.lis
rad_file=$RAD
irr_file=$IRR
cld_file=$CLD

rad_basename=$(basename $rad_file .nc)
irr_basename=$(basename $irr_file .nc)

# Define product file name template
#
lev2_file_fmt=$(mkgranule_name -L 2 -p %s -v $SDPC_PROCESSING_VERSION ${rad_basename}.nc)
lev2_base_fmt=$(basename $lev2_file_fmt .nc)

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

   work_dir_tarfile="${rad_basename}.${work_dir}.tar"
   granule_dir=$(basename $run_dir)

   # We want every product tar file to contain a copy of the archive_subdir file,
   # and we want tar to delete all the files it collects.
   # To avoid deleting the archive_subdir file we create a temporary directory,
   # copy the archive_subdir file over there, and then create the tar file.

   cd $parent_dir
   tmp_dir="$(mktemp -d tmpo3t.XXXXXX)"
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

finish()
{
   tar_product_to_dest_dir "$l2_repro_dir"
}
trap finish EXIT ERR

out_basename=$(printf "$lev2_base_fmt" O3TOT)

etc_dir="$SDPC_PIPE_DIR/etc"

product_dir=.
spectra_dir=.
cloud_dir=.
refdata_dir="$SDPC_REFDATA_DIR/o3_total"

pcf_file="$product_dir/o3_total.pcf"

radiance_file="$rad_file"
irradiance_file="$irr_file"
cloud_file="$cld_file"
product_file="${out_basename}.nc"
template_pcf="${etc_dir}/o3_total/default.pcf.in"

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
 -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
 $template_pcf > $pcf_file

export PGSMSG="${SDPC_ROOT}/msgs"
export PGS_PC_INFO_FILE="$pcf_file"

srun --ntasks=1 --exclusive --output=log_o3_total.txt --job-name=O3TOT \
  L1_o3_total tempo wrt_odl

# SDPTK MET routines litter the directory with temporary files
find . -maxdepth 1 -name "MCFWrite.temp_*" -delete

# Product is finished at this point
insert_fixed_metadata.py $product_file
fix_met_format.py ${product_file}.met
md5sum $product_file > ${product_file}.md5

trap - EXIT
tar_product_to_dest_dir "$l2_out_dir"

case "$rad_basename" in
   *_NRT_* )
     arch_dest="NRT/L2"
     ;;
   * )
     arch_dest="L2"
     ;;
esac

tarfile_path="$l2_out_dir/${rad_basename}.${work_dir}.tar"
archive.sl --delete -a $SDPC_ARCHIVE_DIR -l $arch_dest $tarfile_path
