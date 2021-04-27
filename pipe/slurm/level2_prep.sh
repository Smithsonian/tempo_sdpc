#! /bin/sh
#SBATCH --output=/dev/null

# 0. This script is normally run on a compute node as a batch process
#    to prepare for L2 product generation.
#
# 1. The script is started in a writeable directory, with the
#    command line:
#          $1 = radiance file basename
#          $2 = file defining these variables:
#                rad_file = geolocated radiance file path
#                irr_file = irradiance file path
#                snow_file = path to NSIDC snow and ice cover data file
#
# 2. The first task is to finish processing of the geolocated radiance file
#    by doing the following:
#       - run L1_inr_post
#       - perform wavelength calibration
#       - perform polarization correction
#    The finished L1 radiance file is then stored in the archive.
#
# 3. The second task is generate the L2 cloud product, which is used
#    later on to generate the other Level 2 data products.
#    The finished cloud file is then stored in the archive.
#
# 4. Upon completion of the cloud product, a tar file containing
#    the radiance, irradiance, and cloud products is stored in the
#    L2/incoming directory.
#
# 5. When it's finished, the script cleans up after itself and
#    should leave nothing behind.
#
# On error, a tar file is stored in the L1/repro directory.
#
#---------------------------------------------------------------------

# exit on error
#set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If RADIANCE_POLCORR is not set, define it to be ON.
# To turn off polarization correction, set it to anything else.
: "${RADIANCE_POLCORR:=ON}"

# If SDPC_DIAGNOSTIC_INDEX is not set, define it to be OFF
# To turn on this diagnostic feature, set it to an integer 0 <= n < num_frames_in_granule
: "${SDPC_DIAGNOSTIC_INDEX:=OFF}"

# If USE_FORECAST_MET_DATA is not set, define it to be OFF
# To use forecast data, set it to anything else
: "${USE_FORECAST_MET_DATA:=OFF}"

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_RUN_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-target-file> <file-list-file>"
  exit 1
fi
rad_file="$1"
file_list_file="$2"

# file_list_file should define the following symbols:
#    rad_path
#    irr_file
#    snow_file
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
export PATH="$SDPC_ROOT/bin:$PATH"
etc_dir="$SDPC_ROOT/etc"

l1_out_dir="$SDPC_RUN_DIR/L1/out"
l1_repro_dir="$SDPC_RUN_DIR_MASTER/L1/repro"
l2_incoming="$SDPC_RUN_DIR_MASTER/L2/incoming"
l2_inputs="$SDPC_RUN_DIR_MASTER/L2/inputs"
l2_out_dir="$SDPC_RUN_DIR/L2/out"

# Make a working directory with a local copy of the radiance file.
rad_file_basename=$(basename "$rad_file" .nc)
work_dir="${rad_file_basename}"
/bin/mkdir "$work_dir"
cd $work_dir
/bin/cp "$rad_path" "$rad_file"
/bin/cp "$file_list_file" "${rad_file_basename}.lis"
chmod u+w "$rad_file"

run_dir=$(pwd)
parent_dir=$(dirname "$run_dir")
granule_dir=$(basename "$run_dir")

# Input files:
#
rad_basename=$(basename "$rad_file" .nc)
irr_basename=$(basename "$irr_file" .nc)

# Define template product file name
#
lev1_file_fmt=$(mkgranule_name -L 1 -p %s "${rad_basename}.nc")
lev2_file_fmt=$(mkgranule_name -L 2 -p %s "${rad_basename}.nc")
lev2_base_fmt=$(basename "$lev2_file_fmt" .nc)

# FIXME: For now, both CLDO4 and CLDRR names are needed.
cld_o4_basename=$(printf "$lev2_base_fmt" CLDO4)
cld_rr_basename=$(printf "$lev2_base_fmt" CLDRR)

get_tiepoint_file()
{
   # If the basename has a leading ".", remove it
   rad_path_basename_sans_ext=$(basename $rad_path .nc | sed -e s"/^[.]//")
   rad_path_dir=$(dirname $rad_path)

   # Since the original radiance path is XXXX.Smoothed.nc
   #    we want the tie point file named XXXX.Smoothed.Internal.nc
   tiepoint_path="$rad_path_dir/${rad_path_basename_sans_ext}.Internal.nc"

   if test -f "$tiepoint_path" ; then
      tiepoint_file=$(printf "$lev1_file_fmt" INR)
      /bin/cp $tiepoint_path $tiepoint_file
      chmod u+w $tiepoint_file
   else
      tiepoint_path=""
      tiepoint_file=""
   fi
}

tar_granule_dir_to_dest()
{
   dest_dir=$1
   mkdir -p "$dest_dir"
   cd $parent_dir
   tarfile="${rad_basename}.tar"
   tar c --remove-files -f $dest_dir/.${tarfile} $granule_dir
   /bin/mv $dest_dir/.${tarfile} $dest_dir/${tarfile}
}

finish()
{
   tar_granule_dir_to_dest "$l1_repro_dir"
}
trap finish EXIT ERR

tar_l1_radiance_to_dest()
{
   dest_dir=$1
   mkdir -p "$dest_dir"
   cd $parent_dir

   if test x"$tiepoint_file" != x ; then
      tiepoint_file="$granule_dir/$tiepoint_file"
   fi

   tarfile_rad="${rad_basename}.rad.tar"

   EXTRA_FILES=""
   if test x"$RADIANCE_POLCORR" = x"ON"; then
      EXTRA_FILES="$granule_dir/log_polcorr.txt"
      if test x"$SDPC_DIAGNOSTIC_INDEX" != x"OFF" ; then
         EXTRA_FILES="$EXTRA_FILES $granule_dir/polcorr_*${rad_basename}.nc"
      fi
   fi

   wavecal_logs="$granule_dir/wavecal_logs.tar.gz"
   if test -f $wavecal_logs ; then
      EXTRA_FILES="$EXTRA_FILES $wavecal_logs"
   fi

   tar cf $dest_dir/.${tarfile_rad} \
       $granule_dir/${rad_basename}.nc \
       $granule_dir/archive_subdir \
       $granule_dir/log_inr_post.txt $tiepoint_file \
       $granule_dir/${rad_basename}.nc.met $EXTRA_FILES

   /bin/mv $dest_dir/.${tarfile_rad} $dest_dir/${tarfile_rad}

   archive.sl --clobber --delete -a $SDPC_ARCHIVE_DIR -l L1 $dest_dir/${tarfile_rad}

   /bin/rm -f $granule_dir/log_inr_post.txt $tiepoint_file \
              $granule_dir/${rad_basename}.nc.met $EXTRA_FILES

   # Now that the final L1b radiance file has been archived, we can
   # delete the L1a radiance file that was provided as input to INR:
   /bin/rm -f $SDPC_INR_RUN_DIR/Staging/Granules/${rad_basename}.nc

   # Move INR performance reports to the archive:
   inr_report="$SDPC_INR_RUN_DIR/Output/${rad_basename}.PerformanceReport.nc"
   if test -f "$inr_report" ; then
      scan_dir=$(dirname $SDPC_ARCHIVE_DIR/L1/$granule_subdir)
      inr_dir=$(dirname $scan_dir)/inr
      if ! test -d "$inr_dir" ; then
         mkdir -p $inr_dir
      fi
      /bin/mv $inr_report $inr_dir
   fi
}

. $SDPC_ROOT/bin/wavecal.sh

run_inr_post()
{
   radiance_file=$1

   # INR post
   srun --ntasks=1 --output=log_inr_post.txt \
    L1_inr_post -vv -c ${etc_dir}/l1_inr_post.cfg \
                -s $snow_file $radiance_file

   run_wavecal $radiance_file "0-4"

   # polarization correction
   if test x"$RADIANCE_POLCORR" = x"ON"; then

      # Intentionally avoid TEMPO prefix for diagnostic output files,
      # so those files won't be treated as data products.
      if test x"$SDPC_DIAGNOSTIC_INDEX" != x"OFF" ; then
         copy_level1_frame.py $SDPC_DIAGNOSTIC_INDEX \
               $radiance_file polcorr_before_${rad_basename}.nc
      fi

      # perform polarization correction
      srun --ntasks=1 --output=log_polcorr.txt \
       L1_polcorr -c ${etc_dir}/l1_inr_post.cfg $radiance_file

      if test x"$SDPC_DIAGNOSTIC_INDEX" != x"OFF" ; then
         copy_level1_frame.py $SDPC_DIAGNOSTIC_INDEX \
               $radiance_file polcorr_after_${rad_basename}.nc
      fi
   fi

   # compress before archiving
   num_dsets_gzipped=$(h5stat -d $radiance_file | grep GZIP | cut -d: -f2)
   if test $num_dsets_gzipped -eq 0 ; then
      tmpfile="${radiance_file}.prezip"
      /bin/mv $radiance_file $tmpfile
      srun --ntasks=1 nccopy -s -d 1 $tmpfile $radiance_file
      /bin/rm $tmpfile
   fi

   (tar_l1_radiance_to_dest "$l1_out_dir")
}

tar_l2_cloud_to_dest()
{
   cld_dir=$1
   dest_dir=$2
   mkdir -p "$dest_dir"
   cd $parent_dir
   tarfile_cld="${rad_basename}.cld.tar"
   tar cf $dest_dir/.${tarfile_cld} \
          $granule_dir/archive_subdir $granule_dir/$cld_dir
   /bin/mv $dest_dir/.${tarfile_cld} $dest_dir/${tarfile_cld}

   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $dest_dir/${tarfile_cld}
}

run_cloud_rr()
{
  cld_rr_dir="CLDRR"
  mkdir -p $cld_rr_dir
  cd $cld_rr_dir

  product_dir=.
  spectra_dir=".."
  refdata_dir="$SDPC_RUN_DIR/refdata/cloud"

  radiance_file="${rad_basename}.nc"
  irradiance_file="${irr_basename}.nc"
  product_file="${cld_rr_basename}.nc"
  template_pcf="${etc_dir}/cloud/default.pcf.in"
  pcf_file="$product_dir/cloud.pcf"

  # Edit the PCF file template:
  sed \
      -e s,@refdata_dir@,$refdata_dir,g \
      -e s,@spectra_dir@,$spectra_dir,g \
      -e s,@product_dir@,$product_dir,g \
      -e s,@radiance_file@,$radiance_file,g \
      -e s,@irradiance_file@,$irradiance_file,g \
      -e s,@product_file@,$product_file,g \
      -e s,@etc_dir@,$etc_dir,g \
      -e s,@versionid@,$SDPC_PROCESSING_VERSION,g \
      $template_pcf > $pcf_file

  export PGS_PC_INFO_FILE="$pcf_file"
  export PGSMSG="${SDPC_ROOT}/msgs"

  srun --ntasks=1 --output=log_cloud.txt \
    L1_cloud -tempo -wrt_odl

  # SDPTK MET routines litter the directory with temporary files
  /bin/rm -f MCFWrite.temp_*

  tar_l2_cloud_to_dest "$cld_rr_dir" "$l2_out_dir"

  cd $granule_dir
  /bin/mv ${cld_rr_dir}/${cld_rr_basename}.nc .
  /bin/rm ${cld_rr_dir}/*
  /bin/rmdir ${cld_rr_dir}
}

set_met_file_path()
{
  varname=$1
  met_file_path=$(grep ${varname} ../${rad_basename}.lis | sed -e s,${varname}=,,)
}

derive_o2o2_slant_column()
{
  out_basename=$1
  molecule="O2O2"

  export TG_NO_HE5_OUTPUT=1

  product_dir="."
  spectra_dir=".."
  refsec_dir="$SDPC_RUN_DIR/refdata/trace_gas/refsec"
  refdata_dir="$SDPC_RUN_DIR/refdata"

  # FIXME! reference sector files are in $refsec_dir
  refsec_rad_file="OML1BRUG-o08544.nc"
  refsec_cld_file="OMCLDRR-o08544.nc"

  radiance_file="${rad_basename}.nc"
  irradiance_file="${irr_basename}.nc"
  product_file="${out_basename}.nc"
  template_pcf="${etc_dir}/trace_gas/default.pcf.${molecule}.in"
  template_ctrl="${etc_dir}/trace_gas/control.${molecule}.in"
  control_file="control_${molecule}.txt"
  this_pcf_file="trace_gas.pcf_${molecule}"

  # Select climatology file for the month the data was acquired
  start_month_name=$(level1_info --month ../$radiance_file)
  clim_file="TEMPO_GEOS-Chem_climatology_${start_month_name}_v0p0.he5"

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
   -e s,@product_dir@,$product_dir,g \
   -e s,@etc_dir@,$etc_dir,g \
   -e s,@radiance_file@,$radiance_file,g \
   -e s,@irradiance_file@,$irradiance_file,g \
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
}

run_cloud_o4()
{
  cld_o4_dir="CLDO4"
  mkdir -p $cld_o4_dir
  cd $cld_o4_dir

  # First derive O2O2 slant column using the trace gas code:
  derive_o2o2_slant_column $cld_o4_basename

  # Change labels O2O2 => CLDO4 in both .nc and .met files
  fix_cldo4_metadata.py "${cld_o4_basename}.nc"
  sed -i -e s,TEMPO_O2O2_L2,TEMPO_CLDO4_L2,g "${cld_o4_basename}.nc.met"

  # TBD: Add cloud fields to CLDO4 product file:

  tar_l2_cloud_to_dest "$cld_o4_dir" "$l2_out_dir"

  cd $granule_dir
  /bin/mv ${cld_o4_dir}/${cld_o4_basename}.nc .
  /bin/rm ${cld_o4_dir}/*
  /bin/rmdir ${cld_o4_dir}
}

create_file_listing()
{
  cloud_basename=$1

cat << EOF > pge_input_basenames.lis
RAD=${rad_basename}.nc
IRR=${irr_basename}.nc
CLD=${cloud_basename}.nc
EOF
}

perform_cleanup()
{
   # Delete the original radiance file, and file list file
   /bin/rm "$rad_path" "$file_list_file"

   # We need original radiance path, without the "." prefix on the basename
   # (usually looks like "$some_dir/${prefix}.Smoothed.nc"):
   rad_path_nodot="$rad_path_dir/${rad_path_basename_sans_ext}.nc"

   INR_FILE_TAGS="NavigatedResult Internal Smoothed.Internal"

   for tag in $INR_FILE_TAGS ; do
       inrfile_path=$(echo "$rad_path_nodot" | sed -e "s/Smoothed.nc/${tag}.nc/")
       if test -f "$inrfile_path" ; then
          /bin/rm "$inrfile_path"
       fi
   done
}

notify_granule_ready()
{
  dir="$1"

  tarfile_basename="${rad_basename}.tar"
  this_hostname_sans_domain=$(uname -n | cut -d. -f1)

  # To give notice that a granule is ready for level 1-2 processing,
  # we store path information in an ascii file in $l2_incoming.
  # The file is created as a hidden file, then renamed so that the file
  # appears as an atomic change to the directory file listing.
cat <<EOF > $dir/.$tarfile_basename
tar_host="$this_hostname_sans_domain"
tar_host_file_path="${l1_out_dir}/$tarfile_basename"
granule_arch_dir_path="$SDPC_ARCHIVE_DIR/L2/$granule_subdir"
level2_products="$SDPC_LEVEL2_PRODUCTS"
EOF
  /bin/mv "$dir/.$tarfile_basename" "$dir/$tarfile_basename"
}

if ! test -f "$irr_file" ; then
  echo "ERROR:  irradiance file not found:  $irr_file"
  exit 1
fi

get_tiepoint_file

granule_subdir=$(level1_info --dir ${rad_basename}.nc)
printf "$granule_subdir" > archive_subdir

# We'll be updating the metadata file, so retrieve the pre-INR version from the archive
/bin/cp "$SDPC_ARCHIVE_DIR/L1/$granule_subdir/${rad_basename}.nc.met" .

run_inr_post ${rad_basename}.nc

/bin/cp $irr_file ${irr_basename}.nc
(run_cloud_o4)
(run_cloud_rr)

# FIXME:  For now, we'll delete CLDO4, and use CLDRR for retrievals.
/bin/rm -f "${cld_o4_basename}.nc"
create_file_listing  "$cld_rr_basename"

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

tar_granule_dir_to_dest "$l1_out_dir"

# When level 2 products are to be generated next, put a tar notice file
# in $l2_incoming in preparation for the next processing stage.
# Using a notice file instead of the tar file itself minimizes data
# movement and should improve efficiency.
# When no further processing is intended, put the tar file itself in
# $l2_inputs so that the required input for the next stage is collected
# in one place on the master node.

if ! test x"$SDPC_LEVEL2_PRODUCTS" = x"NONE"; then
   notify_granule_ready "$l2_incoming"
else
   local_tar_file="${l1_out_dir}/${rad_basename}.tar"
   /bin/cp "$local_tar_file" "$l2_inputs"
   /bin/rm "$local_tar_file"
fi

perform_cleanup
