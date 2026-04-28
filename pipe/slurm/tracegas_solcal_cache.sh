#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# 1. Processing ultimately stores all results in a tar file in an
#    appropriate destination directory.
# 2. When processing ends, remove the processing directory.

molecule="$1"
irr_file="$2"

work_dir="$molecule"

l2_out_dir="$SDPC_NODE_DIR/L2/out"
l2_repro_dir="$SDPC_PIPE_DIR/repro/L2"

run_dir=$(pwd)
parent_dir=$(dirname $run_dir)

mkdir -p $work_dir
ln $irr_file $work_dir
cd $work_dir

irr_basename=$(basename $irr_file .nc)

# Define solcal file name template
solcal_file_fmt=$(mkgranule_name -L 2 -p %s -v $SDPC_PROCESSING_VERSION ${irr_basename}.nc)

tar_product_to_dest_dir()
{
   dest_dir=$1

   cd $run_dir

   if ! test -d $work_dir ; then
      return
   fi

   /bin/rm -f $work_dir/${irr_basename}.nc

   work_dir_tarfile="${irr_basename}.${work_dir}.tar"
   granule_dir=$(basename $run_dir)

   # We want every product tar file to contain a copy of the archive_subdir file,
   # and we want tar to delete all the files it collects.
   # To avoid deleting the archive_subdir file we create a temporary directory,
   # copy the archive_subdir file over there, and then create the tar file.

   cd $parent_dir
   tmp_dir="$(mktemp -d tmp${molecule}.XXXXXX)"
   tmp_granule_dir="$tmp_dir/$irr_basename"
   /bin/mkdir -p $tmp_granule_dir
   /bin/cp $granule_dir/archive_subdir $tmp_granule_dir
   /bin/mv $granule_dir/$work_dir $tmp_granule_dir

   cd $tmp_dir

   /bin/mkdir -p $dest_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile \
         $irr_basename/archive_subdir $irr_basename/$work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile

   cd $parent_dir
   /bin/rmdir $tmp_granule_dir $tmp_dir
}

finish()
{
   tar_product_to_dest_dir "$l2_repro_dir"
}
trap finish EXIT ERR

solcal_file=$(printf "$solcal_file_fmt" IRR${molecule})
export TG_NO_HE5_OUTPUT=1

etc_dir="$SDPC_PIPE_DIR/etc"

product_dir=.
spectra_dir=.
cloud_dir=.

pcf_file="$product_dir/trace_gas.pcf"

product_file="notused"
template_pcf="${etc_dir}/trace_gas/default.pcf.${molecule}.in"
template_ctrl="${etc_dir}/trace_gas/control.${molecule}.in"
control_file="control_${molecule}.txt"
this_pcf_file="${pcf_file}_${molecule}"

cld_file=""

# hack for solcal_cache_mode=save (radiance file isn't used, but a valid filename is needed)
rad_file="$irr_file"

radref_basename="notused"
radref_dirname="$spectra_dir"

clim_file=""

met_file1=""
met_dir1=""
met_file2=""
met_dir2=""

# copy the control file template
/bin/cp $template_ctrl $control_file

# Edit the PCF file template:
sed \
 -e s,@control_file@,$control_file,g \
 -e s,@refdata_dir@,$SDPC_REFDATA_DIR,g \
 -e s,@spectra_dir@,$spectra_dir,g \
 -e s,@cloud_dir@,$cloud_dir,g \
 -e s,@product_dir@,$product_dir,g \
 -e s,@solcal_cache_mode@,save,g \
 -e s,@solcal_source@,solar_irradiance,g \
 -e s,@solcal_file@,$solcal_file,g \
 -e s,@solcal_dir@,$spectra_dir,g \
 -e s,@etc_dir@,$etc_dir,g \
 -e s,@radiance_file@,$rad_file,g \
 -e s,@irradiance_file@,$irr_file,g \
 -e s,@cloud_file@,$cld_file,g \
 -e s,@clim_file@,$clim_file,g \
 -e s,@met_dir1@,$met_dir1,g \
 -e s,@met_file1@,$met_file1,g \
 -e s,@met_dir2@,$met_dir2,g \
 -e s,@met_file2@,$met_file2,g \
 -e s,@product_file@,$product_file,g \
 -e s,@radref_basename@,$radref_basename,g \
 -e s,@radref_dirname@,$radref_dirname,g \
 -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
 $template_pcf > $this_pcf_file

export PGS_PC_INFO_FILE="$this_pcf_file"

srun --ntasks=1 --exclusive --quiet --output=log_${molecule}.txt --job-name=${molecule} \
 L1_trace_gas -tempo -wrt_odl

# SDPTK MET routines litter the directory with temporary files
find . -maxdepth 1 -name "MCFWrite.temp_*" -delete

trap - EXIT
tar_product_to_dest_dir "$l2_out_dir"

tarfile_path="$l2_out_dir/${irr_basename}.${work_dir}.tar"
archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
