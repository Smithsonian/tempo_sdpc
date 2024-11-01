#! /bin/sh
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1

# 0. This script is normally run on a compute node as a batch process
#    to prepare for NRT L2 product generation.
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
#       - run L1_inr_post
#
# 3. The second task is to do the O2O2 retrieval as a first step
#    in generating the L2 cloud product, which is used
#    later on to generate the other Level 2 data products.
#    A tar file containing the results is stored on the compute node,
#    and a tar "notice" file is stored in
#              $SDPC_PIPE_DIR/stage/granules/nrt/cldo4_prep
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

# file_list_file should define the following symbols:
#    rad_path
#    irr_file
#    snow_file
#    orig_rad_path
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
etc_dir="$SDPC_PIPE_DIR/etc"

l1_repro_dir="$SDPC_PIPE_DIR/repro/L1"
nrt_incoming="$SDPC_PIPE_DIR/stage/granules/inr_output/nrt/cldo4_prep"
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

# To generate NRT filenames, mkgranule_name needs to find
# the attribute near_real_time=1 in the radiance file header
ncatted -O -h -a near_real_time,global,c,i,1 "$rad_file"

# Define template product file name
#
lev1_file_fmt=$(mkgranule_name -L 1 -p %s -v $SDPC_PROCESSING_VERSION "${rad_basename}.nc")
lev2_file_fmt=$(mkgranule_name -L 2 -p %s -v $SDPC_PROCESSING_VERSION "${rad_basename}.nc")
lev2_base_fmt=$(basename "$lev2_file_fmt" .nc)

cld_o4_basename=$(printf "$lev2_base_fmt" CLDO4)

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
     /bin/rm -f $diagnostic_file
  fi

  # SDPTK MET routines litter the directory with temporary files
  find . -maxdepth 1 -name "MCFWrite.temp_*" -delete
}

tar_l2_cloud_to_dest()
{
   cld_dir=$1
   dest_dir=$2
   mkdir -p "$dest_dir"
   cd $parent_dir
   tarfile_cld="${rad_basename}.cld.tar"
   tar cf $dest_dir/.${tarfile_cld} \
          $granule_dir/archive_subdir $granule_dir/$cld_dir \
          $granule_dir/pge_input_basenames.lis
   /bin/mv $dest_dir/.${tarfile_cld} $dest_dir/${tarfile_cld}

   #archive.sl --delete -a $SDPC_ARCHIVE_DIR -l NRT/L2 $dest_dir/${tarfile_cld}
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

  cd $run_dir
  create_file_listing "$cld_o4_basename"

  tar_l2_cloud_to_dest "$cld_o4_dir" "$l2_out_dir"
}

perform_cleanup()
{
   # Delete the preserved radiance file copy, and file list file
   /bin/rm -f "$rad_path" "$file_list_file"
   /bin/rm -f "$orig_rad_path"
}

notify_o2o2_granule_ready()
{
  dir="$1"

  tarfile_basename="${rad_basename}.cld.tar"
  this_hostname_sans_domain=$(uname -n | cut -d. -f1)

  # To give notice that a granule is ready for level 2 processing,
  # we store path information in an ascii file in $dir.
  # The file is created as a hidden file, then renamed so that the file
  # appears as an atomic change to the directory file listing.
cat <<EOF > $dir/.$tarfile_basename
tar_host="$this_hostname_sans_domain"
tar_host_file_path="${l2_out_dir}/$tarfile_basename"
granule_arch_dir_path="$SDPC_ARCHIVE_DIR/NRT/L2/$granule_subdir"
cld_o4_basename="${cld_o4_basename}.nc"
EOF
  /bin/mv "$dir/.$tarfile_basename" "$dir/$tarfile_basename"
}

if ! test -f "$irr_file" ; then
  echo "ERROR:  irradiance file not found:  $irr_file"
  exit 1
fi

granule_subdir=$(level1_info --dir ${rad_basename}.nc)
printf "$granule_subdir" > archive_subdir

run_inr_post ${rad_basename}.nc

/bin/cp $irr_file ${irr_basename}.nc
(run_cloud_o4)

cd $parent_dir
/bin/rm -rf $run_dir

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# To facilitate generation of level 2 products, put a tar notice file
# in $nrt_incoming. Using a notice file instead of the tar
# file itself minimizes data movement and should improve efficiency.

notify_o2o2_granule_ready "$nrt_incoming"

perform_cleanup
