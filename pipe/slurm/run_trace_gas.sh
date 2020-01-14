#! /bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If USE_FORECAST_MET_DATA is not set, define it to be OFF
# To use forecast data, set it to anything else
: "${USE_FORECAST_MET_DATA:=OFF}"

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
. ./pge_input_basenames.lis
rad_file=$RAD
irr_file=$IRR
cld_file=$CLD

rad_basename=$(basename $rad_file .nc)
irr_basename=$(basename $irr_file .nc)

# Define product file name template
#
lev2_file_fmt=$(mkgranule_name -L 2 -p %s ${rad_basename}.nc)
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

product_dir=.
spectra_dir=.
cloud_dir=.
refsec_dir="$SDPC_RUN_DIR/refdata/trace_gas/refsec"
refdata_dir="$SDPC_RUN_DIR/refdata"

pcf_file="$product_dir/trace_gas.pcf"

# FIXME! reference sector files are in $refsec_dir
refsec_rad_file="OML1BRUG-o08544.nc"
refsec_cld_file="OMCLDRR-o08544.nc"

product_file="${out_basename}.nc"
template_pcf="${etc_dir}/trace_gas/synth.pcf.${molecule}.in"
template_ctrl="${etc_dir}/trace_gas/synth.control.${molecule}.in"
control_file="control_${molecule}.txt"
this_pcf_file="${pcf_file}_${molecule}"

# Select climatology file for the month the data was acquired
declare -a month_names=(INVALID JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)
start_month=$(grep tstart_month granule_ident.csv | cut -f2 -d,)
start_month_name=${month_names[${start_month}]}
clim_file="TEMPO_GEOS-Chem_climatology_${start_month_name}_v0p0.he5"

set_met_file_path()
{
  varname=$1
  met_file_path=$(grep ${varname} ${rad_basename}.lis | sed -e s,${varname}=,,)
}

if test x"$USE_FORECAST_MET_DATA" = x"OFF" ; then
      set_met_file_path "met_file_path_synth"
      met_file1=$(basename $met_file_path)
      met_dir1=$(dirname $met_file_path)

      met_file2=""
      met_dir2=""
else
      set_met_file_path "met_file_path_hires"
      met_file1=$(basename $met_file_path)
      met_dir1=$(dirname $met_file_path)

      set_met_file_path "met_file_path_lores"
      met_file2=$(basename $met_file_path)
      met_dir2=$(dirname $met_file_path)
fi

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
 -e s,@refsec_rad_file@,$refsec_rad_file,g \
 -e s,@refsec_cld_file@,$refsec_cld_file,g \
 -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
 $template_pcf > $this_pcf_file

export PGS_PC_INFO_FILE="$this_pcf_file"

srun --ntasks=1 --output=log_${molecule}.txt \
 L1_trace_gas -tempo -wrt_odl

# SDPTK MET routines litter the directory with temporary files
/bin/rm -f MCFWrite.temp_*

trap - EXIT
tar_product_to_dest_dir "$l2_out_dir"

tarfile_path="$l2_out_dir/${rad_basename}.${work_dir}.tar"
archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
