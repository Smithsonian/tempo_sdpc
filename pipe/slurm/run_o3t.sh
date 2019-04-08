#! /bin/sh
#SBATCH --cpus-per-task=1
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

work_dir="o3t"

l2_out_dir="$SDPC_RUN_DIR/L2/out"
l2_repro_dir="$SDPC_RUN_DIR/L2/repro"

run_dir=$(pwd)
parent_dir=$(dirname $run_dir)
cd $work_dir

# get input file names
. ./files.lis
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

   /bin/rm $work_dir/${rad_basename}.nc
   /bin/rm $work_dir/${irr_basename}.nc
   /bin/rm $work_dir/$cld_file

   work_dir_tarfile="${rad_basename}.${work_dir}.tar"
   granule_dir=$(basename $run_dir)

   /bin/mkdir -p $dest_dir
   /bin/mv $work_dir/granule_ident.csv .
   cd $parent_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile \
         $granule_dir/granule_ident.csv $granule_dir/$work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile
}

finish()
{
   tar_product_to_dest_dir "$l2_repro_dir"
}
trap finish EXIT ERR

out_basename=$(printf "$lev2_base_fmt" o3t)

etc_dir="$SDPC_ROOT/etc"

config_file="$etc_dir/o3_total/o3_total.rc"
. $config_file

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
 $template_pcf > $pcf_file

export PGSMSG="${SDPC_ROOT}/msgs"
export PGS_PC_INFO_FILE="$pcf_file"

srun --ntasks=1 --output=log_o3_total.txt \
  L1_o3_total tempo wrt_odl

# SDPTK MET routines litter the directory with temporary files
/bin/rm -f MCFWrite.temp_*

trap - EXIT
tar_product_to_dest_dir "$l2_out_dir"

tarfile_path="$l2_out_dir/${rad_basename}.${work_dir}.tar"
archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
