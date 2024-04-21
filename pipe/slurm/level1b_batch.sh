#! /bin/sh
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1

# 0. This script is normally run on a compute node as a batch process
#    to prepare for L2 product generation.
#
# 1. The script is started in a writeable directory, with the
#    command line:
#          $1 = radiance file basename
#          $2 = file defining these variables:
#                rad_file = geolocated radiance file path (or RADT_L1)
#                irr_file = irradiance file path
#                snow_file = path to NSIDC snow and ice cover data file
#
# 2. The first task is to finish processing of the geolocated radiance file
#    by doing the following:
#       - run L1_inr_post  [only step done for RADT]
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
#    stage/granules/level2_input directory.
#
# 5. When it's finished, the script cleans up after itself and
#    should leave nothing behind.
#
# On error, a tar file is stored in the repro/L1 directory.
#
#---------------------------------------------------------------------

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If SDPC_RADIANCE_POSTINR is not set, define it to be ON (non-zero).
# To turn off post-INR radiance processing, set it to zero.
: "${SDPC_RADIANCE_POSTINR:=1}"

# If SDPC_RADIANCE_WAVECAL is not set, define it to be ON (non-zero).
# To turn off radiance wavelength calibration, set it to zero.
: "${SDPC_RADIANCE_WAVECAL:=1}"
: "${SDPC_RADIANCE_WAVECAL_NTASKS:=2}"

# If SDPC_RADIANCE_POLCORR is not set, define it to be ON (non-zero).
# To turn off polarization correction, set it to zero
: "${SDPC_RADIANCE_POLCORR:=1}"

# If SDPC_DIAGNOSTIC_INDEX is not set, define it to be OFF
# To turn on this diagnostic feature, set it to an integer 0 <= n < num_frames_in_granule
: "${SDPC_DIAGNOSTIC_INDEX:=OFF}"

# If USE_SYNTHETIC_MET_DATA is not set, define it to be OFF
# To use synthetic met data, set it to the relevant file path
: "${USE_SYNTHETIC_MET_DATA:=OFF}"

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_NODE_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-target-file> <file-list-file>"
  exit 1
fi
rad_file="$1"
file_list_file="$2"

# Set a flag to indicate this is RADT_L1
is_radt_l1=0
case "$rad_file" in
 *_RADT_* )
     is_radt_l1=1
     ;;

 * )
     ;;
esac

# file_list_file should define the following symbols:
#    rad_path
#    irr_file
#    snow_file
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
etc_dir="$SDPC_PIPE_DIR/etc"

l1_out_dir="$SDPC_NODE_DIR/L1/out"
l1_repro_dir="$SDPC_PIPE_DIR/repro/L1"
l2_incoming="$SDPC_PIPE_DIR/stage/granules/level2_input"
l2_out_dir="$SDPC_NODE_DIR/L2/out"

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
lev1_file_fmt=$(mkgranule_name -L 1 -p %s -v $SDPC_PROCESSING_VERSION "${rad_basename}.nc")
lev2_file_fmt=$(mkgranule_name -L 2 -p %s -v $SDPC_PROCESSING_VERSION "${rad_basename}.nc")
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
   if test -d $granule_dir ; then
      tarfile="${rad_basename}.tar"
      tar c --remove-files -f $dest_dir/.${tarfile} $granule_dir
      /bin/mv $dest_dir/.${tarfile} $dest_dir/${tarfile}
   fi
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
   if test $SDPC_RADIANCE_POLCORR -ne 0 ; then
      EXTRA_FILES="$granule_dir/log_polcorr.txt"
      if test x"$SDPC_DIAGNOSTIC_INDEX" != x"OFF" ; then
         EXTRA_FILES="$EXTRA_FILES $granule_dir/polcorr_*${rad_basename}.nc"
      fi
   fi

   wavecal_files="wavecal_logs.tar.gz wavecal_joblog.out log_wavecal_merge.txt"
   for f in $wavecal_files ; do
       f_path="$granule_dir/$f"
       if test -f $f_path ; then
          EXTRA_FILES="$EXTRA_FILES $f_path"
       fi
   done

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
   # delete the L1a radiance file that was provided as input to INR,
   # along with any earlier telemetry-only radiance files.
   inr_input_cache="$SDPC_PIPE_DIR/inr/Staging/Granules"
   level1a_granule_path="${inr_input_cache}/${rad_basename}.nc"
   #*******************************************************************
   # As of INRSW R2.3.5, the INR SW supports a warm restart to occur
   # at the start of each day, and following any spacecraft maneuvers
   # that interrupt radiance scanning.
   # This feature requires that the INR SW assume responsibility for
   # deleting files from the INR input cache.  Henceforth, the SDPC pipeline
   # must not delete radiance granules from the INR input cache.
   # The telemetry-only granules are probably obsolete now as well,
   # but we'll leave support for that in place just in case.
   #*******************************************************************
   #if test -f "$level1a_granule_path" ; then
   #   radiance_telem_only.py --delete --before "$level1a_granule_path" "$inr_input_cache"
   #   /bin/rm -f "$level1a_granule_path"
   #fi

   # Move INR performance reports to the archive:
   inr_report="$SDPC_PIPE_DIR/inr/Output/${rad_basename}.PerformanceReport.nc"
   if test -f "$inr_report" ; then
      scan_dir=$(dirname $SDPC_ARCHIVE_DIR/L1/$granule_subdir)
      inr_dir=$(dirname $scan_dir)/inr
      if ! test -d "$inr_dir" ; then
         mkdir -p $inr_dir
      fi
      /bin/mv $inr_report $inr_dir
   fi
}

run_inr_post()
{
   radiance_file=$1

   # INR post
   srun --ntasks=1 --output=log_inr_post.txt \
    L1_inr_post -vv -c ${etc_dir}/l1_inr_post.cfg \
                -s $snow_file $radiance_file

   domain_file="${etc_dir}/asdc_collection_domain.csv"
   if test -f $domain_file ; then
      inr_ok=$(polygon_inside_domain.py --domain $domain_file $radiance_file)
      if test x"$inr_ok" != xyes ; then
         echo "ERROR: bounding polygon outside collection domain: $radiance_file"
         exit 1
      fi
   fi

   # This batch script must specify ntasks >= 2x
   # the number of wavecal tasks specified here
   if test $SDPC_RADIANCE_WAVECAL -ne 0 -a $is_radt_l1 -eq 0 ; then
      wavecal.sh $radiance_file $SDPC_RADIANCE_WAVECAL_NTASKS
   fi

   # polarization correction
   if test $SDPC_RADIANCE_POLCORR -ne 0 -a $is_radt_l1 -eq 0 ; then

      # Intentionally avoid TEMPO prefix for diagnostic output files,
      # so those files won't be treated as data products.
      if test x"$SDPC_DIAGNOSTIC_INDEX" != x"OFF" ; then
         copy_level1_frame.py $SDPC_DIAGNOSTIC_INDEX \
               $radiance_file polcorr_before_${rad_basename}.nc
      fi

      # perform polarization correction
      srun --ntasks=1 --output=log_polcorr.txt \
       L1_polcorr -c ${etc_dir}/l1_polcorr.cfg $radiance_file

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

   # If the $input_files entry in the RAD_L1 .met file hasn't
   # been expanded yet, expand it here.
   fixup_radl1_inputfiles.py $radiance_file

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
  refdata_dir="$SDPC_REFDATA_DIR/cloud"

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
  find . -maxdepth 1 -name "MCFWrite.temp_*" -delete

  tar_l2_cloud_to_dest "$cld_rr_dir" "$l2_out_dir"

  cd $granule_dir
  /bin/mv ${cld_rr_dir}/${cld_rr_basename}.nc .
  /bin/rm ${cld_rr_dir}/*
  /bin/rmdir ${cld_rr_dir}
}

derive_o2o2_slant_column()
{
  out_basename=$1
  molecule="O2O2"

  export TG_NO_HE5_OUTPUT=1

  product_dir="."
  spectra_dir=".."
  refsec_dir="$SDPC_REFDATA_DIR/trace_gas/refsec"

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
   -e s,@refdata_dir@,"$SDPC_REFDATA_DIR",g \
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

  diagnostic_file="${out_basename}_diag.nc"
  if test -f $diagnostic_file ; then
     tg_diag_filter.py --output diaglog_${molecule}.nc log_${molecule}.txt $diagnostic_file
  fi

  # SDPTK MET routines litter the directory with temporary files
  find . -maxdepth 1 -name "MCFWrite.temp_*" -delete
}

derive_cloud_o4_params()
{
   product_file="$1"

   radiance_file="${rad_basename}.nc"
   irradiance_file="${irr_basename}.nc"

   refdata_dir="$SDPC_REFDATA_DIR/cloud_o4"
   template_ctrl="${etc_dir}/cloud_o4/control.txt.in"
   ctrl_file="cloud_o4_control.txt"

   # edit the control file template
   sed -e "s,@cldo4_file@,$product_file," \
       -e "s,@radiance_file@,$radiance_file," \
       -e "s,@irradiance_file@,$irradiance_file," \
       -e "s,@product_file@,$product_file," \
       -e "s,@refdata_dir@,$refdata_dir," \
       $template_ctrl > $ctrl_file

   srun --ntasks=1 --output=log_cloud_o4.txt \
        L1_cloud_o4 $ctrl_file

   # remove variables:
   #        - main_data_quality_flag
   #        - surface_pressure
   VARS_TO_REMOVE="main_data_quality_flag surface_pressure"
   for var in $VARS_TO_REMOVE ; do
       ncks -x -v $var $product_file ${product_file}.tmp
       /bin/mv ${product_file}.tmp $product_file
   done
   # rename SurfacePressure -> surface_pressure
   ncrename -v SurfacePressure,surface_pressure $product_file
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

  # From O2O2 slant column, derive cloud parameters
  derive_cloud_o4_params "${cld_o4_basename}.nc"

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
   /bin/rm -f "$rad_path" "$file_list_file"

   # We need original radiance path, without the "." prefix on the basename
   # (usually looks like "$some_dir/${prefix}.Smoothed.nc"):
   rad_path_nodot="$rad_path_dir/${rad_path_basename_sans_ext}.nc"

   INR_FILE_TAGS="NavigatedResult Internal Smoothed.Internal"

   for tag in $INR_FILE_TAGS ; do
       inrfile_path=$(echo "$rad_path_nodot" | sed -e "s/Smoothed.nc/${tag}.nc/")
       if test -f "$inrfile_path" ; then
          /bin/rm -f "$inrfile_path"
       fi
   done
}

notify_granule_ready()
{
  dir="$1"

  tarfile_basename="${rad_basename}.tar"
  this_hostname_sans_domain=$(uname -n | cut -d. -f1)

  # To give notice that a granule is ready for level 1-2 processing,
  # we store path information in an ascii file in $dir.
  # The file is created as a hidden file, then renamed so that the file
  # appears as an atomic change to the directory file listing.
cat <<EOF > $dir/.$tarfile_basename
tar_host="$this_hostname_sans_domain"
tar_host_file_path="${l1_out_dir}/$tarfile_basename"
granule_arch_dir_path="$SDPC_ARCHIVE_DIR/L2/$granule_subdir"
rad_filename="${rad_basename}.nc"
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

if test $SDPC_RADIANCE_POSTINR -ne 0 ; then
   # We'll be updating the metadata file, so retrieve the pre-INR version from the archive
   /bin/cp "$SDPC_ARCHIVE_DIR/L1/$granule_subdir/${rad_basename}.nc.met" .
   run_inr_post ${rad_basename}.nc
fi

if test $is_radt_l1 -eq 0 ; then
   /bin/cp $irr_file ${irr_basename}.nc
   (run_cloud_o4)
   create_file_listing  "$cld_o4_basename"
fi

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# To facilitate generation of level 2 products, put a tar notice file
# in $l2_incoming. Using a notice file instead of the tar
# file itself minimizes data movement and should improve efficiency.

if test $is_radt_l1 -eq 0 ; then
   tar_granule_dir_to_dest "$l1_out_dir"
   notify_granule_ready "$l2_incoming"
fi

perform_cleanup
