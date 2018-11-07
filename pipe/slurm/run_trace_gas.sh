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

molecule="$1"
work_dir="$molecule"

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

out_basename=$(printf "$lev2_base_fmt" ${molecule})
export TG_NO_HE5_OUTPUT=1

etc_dir="$SDPC_ROOT/etc"

config_file="$etc_dir/trace_gas/trace_gas.rc"
. $config_file

radiance_file="$rad_file"
irradiance_file="$irr_file"
cloud_file="$cld_file"
product_file="${out_basename}.he5"
template_pcf="${etc_dir}/trace_gas/synth.pcf.${molecule}.in"
template_ctrl="${etc_dir}/trace_gas/synth.control.${molecule}.in"
control_file="control_${molecule}.txt"
this_pcf_file="${pcf_file}_${molecule}"

# Select climatology file for the month the data was acquired
declare -a month_names=(INVALID JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)
start_month=$(grep tstart_month granule_ident.csv | cut -f2 -d,)
start_month_name=${month_names[${start_month}]}
clim_file="TEMPO_GEOS-Chem_climatology_${start_month_name}_v0p0.he5"

# read file-list file to obtain definition for met_file_path variable
met_file_path=$(grep met_file_path ${rad_basename}.lis | sed -e s,met_file_path=,,)

# set meteorological data file:
met_dir=$(dirname $met_file_path)
met_file=$(basename $met_file_path)

# copy the control file template
/bin/cp $template_ctrl $control_file

# Edit the PCF file template:
sed \
 -e s,@control_file@,$control_file,g \
 -e s,@refdata_dir@,$refdata_dir,g \
 -e s,@refsec_dir@,$refsec_dir,g \
 -e s,@spectra_dir@,$spectra_dir,g \
 -e s,@cloud_dir@,$cloud_dir,g \
 -e s,@product_dir@,$product_dir,g \
 -e s,@radiance_file@,$radiance_file,g \
 -e s,@irradiance_file@,$irradiance_file,g \
 -e s,@cloud_file@,$cloud_file,g \
 -e s,@clim_file@,$clim_file,g \
 -e s,@met_dir@,$met_dir,g \
 -e s,@met_file@,$met_file,g \
 -e s,@product_file@,$product_file,g \
 -e s,@refsec_rad_file@,$refsec_rad_file,g \
 -e s,@refsec_cld_file@,$refsec_cld_file,g \
 $template_pcf > $this_pcf_file

export PGS_PC_INFO_FILE="$this_pcf_file"

srun --ntasks=1 --output=log_${molecule}.txt \
 L1_trace_gas -tempo

trap - EXIT
tar_product_to_dest_dir "$l2_out_dir"

tarfile_path="$l2_out_dir/${rad_basename}.${work_dir}.tar"
archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
