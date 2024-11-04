#! /bin/sh
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If USE_SYNTHETIC_MET_DATA is not set, define it to be OFF
# To use synthetic met data, set it to the relevant file path
: "${USE_SYNTHETIC_MET_DATA:=OFF}"

# 1. Processing ultimately stores all results in a tar file in an
#    appropriate destination directory.
# 2. When processing ends, remove the processing directory.

radref_file=""
destripe_file=""

molecule="$1"
shift
if test $# -ge 1 ; then
   radref_file="$1"
   shift
fi
if test $# -ge 1 ; then
   destripe_file="$1"
   shift
fi
work_dir="$molecule"

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

   /bin/rm -f $work_dir/${rad_basename}.nc
   /bin/rm -f $work_dir/${irr_basename}.nc
   /bin/rm -f $work_dir/$cld_file

   work_dir_tarfile="${rad_basename}.${work_dir}.tar"
   granule_dir=$(basename $run_dir)

   # We want every product tar file to contain a copy of the archive_subdir file,
   # and we want tar to delete all the files it collects.
   # To avoid deleting the archive_subdir file we create a temporary directory,
   # copy the archive_subdir file over there, and then create the tar file.

   cd $parent_dir
   tmp_dir="$(mktemp -d tmp${molecule}.XXXXXX)"
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

out_basename=$(printf "$lev2_base_fmt" ${molecule})
export TG_NO_HE5_OUTPUT=1

etc_dir="$SDPC_PIPE_DIR/etc"

product_dir=.
spectra_dir=.
cloud_dir=.
refsec_dir="$SDPC_REFDATA_DIR/trace_gas/refsec"

pcf_file="$product_dir/trace_gas.pcf"

# FIXME! reference sector files are in $refsec_dir
refsec_rad_file="OML1BRUG-o08544.nc"
refsec_cld_file="OMCLDRR-o08544.nc"

product_file="${out_basename}.nc"
template_pcf="${etc_dir}/trace_gas/default.pcf.${molecule}.in"
template_ctrl="${etc_dir}/trace_gas/control.${molecule}.in"
control_file="control_${molecule}.txt"
this_pcf_file="${pcf_file}_${molecule}"

radref_basename=""
radref_dirname=""
case "$molecule" in
  HCHO )
     if test -n "$radref_file" ; then
        radref_basename="$(basename $radref_file)"
        radref_dirname="$(dirname $radref_file)"
     fi
    ;;
  *)
    ;;
esac

# Select climatology file for the month the data was acquired
start_month_name=$(level1_info --month ${rad_basename}.nc)
clim_file="TEMPO_GEOS-Chem_climatology_${start_month_name}_v0p0.he5"

met_file1=""
met_dir1=""
met_file2=""
met_dir2=""
if ! test x"$USE_SYNTHETIC_MET_DATA" = x"OFF" ; then
      met_file1=$(basename $USE_SYNTHETIC_MET_DATA)
      met_dir1=$(dirname $USE_SYNTHETIC_MET_DATA)
fi

# copy the control file template
/bin/cp $template_ctrl $control_file

# Edit the PCF file template:
sed \
 -e s,@control_file@,$control_file,g \
 -e s,@refdata_dir@,$SDPC_REFDATA_DIR,g \
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
 -e s,@radref_basename@,$radref_basename,g \
 -e s,@radref_dirname@,$radref_dirname,g \
 -e s,@refsec_rad_file@,$refsec_rad_file,g \
 -e s,@refsec_cld_file@,$refsec_cld_file,g \
 -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
 $template_pcf > $this_pcf_file

export PGS_PC_INFO_FILE="$this_pcf_file"

srun --ntasks=1 --exclusive --output=log_${molecule}.txt --job-name=${molecule} \
 L1_trace_gas -tempo -wrt_odl

if test -n "$destripe_file" ; then
   case "$molecule" in
         * )
           ;;
      HCHO )
           destripe.py --corrfile $destripe_file $product_file > log_destripe.txt 2>&1
           ;;
   esac
fi

# SDPTK MET routines litter the directory with temporary files
find . -maxdepth 1 -name "MCFWrite.temp_*" -delete

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
