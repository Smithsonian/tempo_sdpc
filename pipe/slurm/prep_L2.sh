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
# On error, a tar file is stored in the L2/repro directory.
#
#---------------------------------------------------------------------

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# If RADIANCE_POLCORR is not set, define it to be ON.
# To turn off polarization correction, set it to anything else.
: "${RADIANCE_POLCORR:=ON}"

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
l2_incoming="$SDPC_RUN_DIR/L2/incoming"
l2_out_dir="$SDPC_RUN_DIR/L2/out"
l2_repro_dir="$SDPC_RUN_DIR/L2/repro"

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
cld_basename=$(printf "$lev2_base_fmt" CLDRR)
cld_dir="CLDRR"

get_tiepoint_file()
{
   # If the basename has a leading ".", remove it
   rad_path_basename_sans_ext=$(basename $rad_path .nc | sed -e s"/^[.]//")
   rad_path_dir=$(dirname $rad_path)

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
   tar_granule_dir_to_dest "$l2_repro_dir"
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
   fi

   tar cf $dest_dir/.${tarfile_rad} \
       $granule_dir/${rad_basename}.nc \
       $granule_dir/granule_ident.csv \
       $granule_dir/log_inr_post.txt $tiepoint_file \
       $granule_dir/${rad_basename}.nc.met $EXTRA_FILES

   /bin/mv $dest_dir/.${tarfile_rad} $dest_dir/${tarfile_rad}

   archive.sl --clobber --delete -a $SDPC_ARCHIVE_DIR -l L1 $dest_dir/${tarfile_rad}

   /bin/rm -f $granule_dir/log_inr_post.txt $tiepoint_file \
              $granule_dir/${rad_basename}.nc.met $EXTRA_FILES
}

. $SDPC_ROOT/bin/run_wavecal.sh

run_inr_post()
{
   radiance_file=$1

   # INR post
   srun --ntasks=1 --output=log_inr_post.txt \
    L1_inr_post -c ${etc_dir}/l1_inr_post.cfg \
                -s $snow_file $radiance_file

   run_wavecal $radiance_file "0-4"

   # polarization correction
   if test x"$RADIANCE_POLCORR" = x"ON"; then
      srun --ntasks=1 --output=log_polcorr.txt \
       L1_polcorr -c ${etc_dir}/l1_inr_post.cfg $radiance_file
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
   dest_dir=$1
   mkdir -p "$dest_dir"
   cd $parent_dir
   tarfile_cld="${rad_basename}.cld.tar"
   tar cf $dest_dir/.${tarfile_cld} \
          $granule_dir/granule_ident.csv $granule_dir/$cld_dir
   /bin/mv $dest_dir/.${tarfile_cld} $dest_dir/${tarfile_cld}

   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $dest_dir/${tarfile_cld}
}

run_cloud()
{
  mkdir -p $cld_dir
  cd $cld_dir

  product_dir=.
  spectra_dir=.
  refdata_dir="$SDPC_RUN_DIR/refdata/cloud"

  radiance_file="../${rad_basename}.nc"
  irradiance_file="../${irr_basename}.nc"
  product_file="${cld_basename}.nc"
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

  tar_l2_cloud_to_dest "$l2_out_dir"

  cd $granule_dir
  /bin/mv ${cld_dir}/${cld_basename}.nc .
  /bin/rm ${cld_dir}/*
  /bin/rmdir ${cld_dir}
}

create_file_listing()
{
cat << EOF > pge_input_basenames.lis
RAD=${rad_basename}.nc
IRR=${irr_basename}.nc
CLD=${cld_basename}.nc
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

get_metadata_file()
{
   granule_ident=$(cat granule_ident.csv | tr ',' '=')
   eval "$granule_ident"

   arch_type="$SDPC_ARCHIVE_DIR/L1/$processing_version/$product_type"
   granule_arch_dir_path="${arch_type}/${sat_local_day_start}/${scan_num}/${granule_num}"

   /bin/cp $granule_arch_dir_path/${rad_basename}.nc.met .
}

if ! test -f "$irr_file" ; then
  echo "ERROR:  irradiance file not found:  $irr_file"
  exit 1
fi

get_tiepoint_file

mkgranule_ident -o granule_ident.csv ${rad_basename}.nc

# We'll be updating the metadata file, so retrieve the pre-INR version from the archive
(get_metadata_file)

run_inr_post ${rad_basename}.nc

/bin/cp $irr_file ${irr_basename}.nc
(run_cloud)

create_file_listing

trap - EXIT
tar_granule_dir_to_dest "$l2_incoming"

perform_cleanup
