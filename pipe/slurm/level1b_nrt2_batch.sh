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
#                rad_path = geolocated radiance file path
#                irr_file = irradiance file path
#                snow_file = path to NSIDC snow and ice cover data file
#    NOTE: The basename ($1) is the desired final basename and need not
#          be the same as $(basename $rad_file)
#
# 2. The first task is to finish processing of the geolocated radiance file
#    by doing the following:
#       - run L1_inr_post
#
# 3. The second task is to finish generating the L2 cloud product,
#    which is used later on to generate the other Level 2 data products.
#    A tar file containing the results is stored on the compute node,
#    and a tar "notice" file is stored in
#              $SDPC_PIPE_DIR/stage/granules/level2_input_nrt
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
#    solcal_file_list
#    destripe_file_list
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
etc_dir="$SDPC_PIPE_DIR/etc"

l1_out_dir="$SDPC_NODE_DIR/L1/out"
l1_repro_dir="$SDPC_PIPE_DIR/repro/L1"
l2_incoming_nrt="$SDPC_PIPE_DIR/stage/granules/level2_input_nrt"
l2_out_dir="$SDPC_NODE_DIR/L2/out"

# Construct the path to the original RAD_L1a .met file because we'll need that later.
orig_granule_subdir=$(level1_info --dir $rad_path)
orig_rad_basename="$(basename $rad_path .Smoothed.nc | sed -e 's/^[.]//' -e 's/_NRT//')"
orig_metfile_path="$SDPC_ARCHIVE_DIR/L1/${orig_granule_subdir}/${orig_rad_basename}.nc.met"

# Input files:
#
rad_basename=$(basename "$rad_file" .nc)
irr_basename=$(basename "$irr_file" .nc)

# Read notice file to locate result from level1b_nrt1.
# This defines the following variables:
# tar_host
# tar_host_file_path
# granule_arch_dir_path
# cld_o4_basename
#
cldo4_input_dir="$SDPC_PIPE_DIR/stage/granules/cldo4_input_nrt"
. $cldo4_input_dir/${rad_basename}.cld.tar

cld_o4_basename=$(basename "$cld_o4_basename" .nc)

# Retrieve the tar file and unpack it:
tar_file_basename=$(basename $tar_host_file_path)
this_host=$(uname -n | cut -d. -f1)
if test x"$tar_host" != x"$this_host" ; then
   scp -o StrictHostKeyChecking=no $tar_host:$tar_host_file_path .
   tar xf $tar_file_basename
   tar_file_dir="$(tar tf $tar_file_basename | head -n1 | cut -d/ -f1)"
   /bin/rm $tar_file_basename
else
   tar xf $tar_host_file_path
   tar_file_dir="$(tar tf $tar_host_file_path | head -n1 | cut -d/ -f1)"
fi

# tar directory will be the working directory

work_dir="${tar_file_dir}"
cd $work_dir
/bin/cp "$rad_path" "$rad_file"
/bin/cp "$file_list_file" "${rad_basename}.lis"
if test -s "$solcal_file_list" ; then
   # We want this list copied into the tar file
   # that provides inputs for all L2 product generation
   /bin/cp "$solcal_file_list" .
fi
if test -s "$destripe_file_list" ; then
   # We want this list copied into the tar file
   # that provides inputs for all L2 product generation
   /bin/cp "$destripe_file_list" .
fi
chmod u+w "$rad_file"

run_dir=$(pwd)
parent_dir=$(dirname "$run_dir")
granule_dir=$(basename "$run_dir")

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

   tarfile_rad="${rad_basename}.rad.tar"

   tar cf $dest_dir/.${tarfile_rad} \
       $granule_dir/${rad_basename}.nc \
       $granule_dir/archive_subdir \
       $granule_dir/log_inr_post.txt \
       $granule_dir/${rad_basename}.nc.met

   /bin/mv $dest_dir/.${tarfile_rad} $dest_dir/${tarfile_rad}

   archive.sl --clobber --delete -a $SDPC_ARCHIVE_DIR -l NRT/L1 $dest_dir/${tarfile_rad}

   /bin/rm -f $granule_dir/log_inr_post.txt \
              $granule_dir/${rad_basename}.nc.met
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
   fix_radl1_metadata.py $radiance_file

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

   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l NRT/L2 $dest_dir/${tarfile_cld}
}

derive_cloud_o4_params()
{
   product_file="$1"

   radiance_file="${rad_basename}.nc"
   irradiance_file="${irr_basename}.nc"

   refdata_dir="$SDPC_REFDATA_DIR/cloud_o4"
   template_ctrl="${etc_dir}/cloud_o4/control.txt.in"
   ctrl_file="cloud_o4_control.txt"

   out_basename="$(basename $product_file .nc)"

   destripe_file="nonexistent"
   if test -s "$destripe_file_list" ; then
      destripe_file_from_lookup="$(grep DSTRCLDO4 $destripe_file_list || true)"
      if test -f "$destripe_file_from_lookup" ; then
         destripe_file="$destripe_file_from_lookup"
      fi
   fi
   # Destripe if possible, but nonexistent destripe file is ok here:
   tempo_destripe_regular.py --mode apply --desfnm "$destripe_file" --l2fnm "$product_file" > log_destripe.txt 2>&1

   # edit the control file template
   sed -e "s,@cldo4_file@,$product_file," \
       -e "s,@radiance_file@,$radiance_file," \
       -e "s,@irradiance_file@,$irradiance_file," \
       -e "s,@product_file@,$product_file," \
       -e "s,@refdata_dir@,$refdata_dir," \
       -e "s,@apply_destripe@,1," \
       -e "s,@apply_solshift@,0," \
       -e "s,@apply_radshift@,0," \
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
  # CLDO4 directory already contains the output from the
  # O2O2 slant column retrieval.
  cld_o4_dir="CLDO4"
  cd $cld_o4_dir

  # Copy the 2nd pass geolocation variables from the radiance
  # file to the cloud data product
  # (L1_cloud_o4 will recompute relative_azimuth_angle)
  update_cldo4_geovars.py "../${rad_basename}.nc" "${cld_o4_basename}.nc"

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
   # Delete the preserved radiance file copy, and file list file
   /bin/rm -f "$rad_path" "$file_list_file" "$solcal_file_list" "$destripe_file_list"
   # Delete this tar notice file and the corresponding tar file
   /bin/rm "$cldo4_input_dir/${rad_basename}.cld.tar"
   if test x"$tar_host" != x"$this_host" ; then
      ssh -o ForwardX11=no "$tar_host" /bin/rm -f "$tar_host_file_path"
   else
      /bin/rm -f "$tar_host_file_path"
   fi
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

# Define archive_subdir before we need to archive anything
granule_subdir=$(level1_info --dir ${rad_basename}.nc)
printf "$granule_subdir" > archive_subdir

# To generate NRT filenames downstream, mkgranule_name needs to find
# the attribute near_real_time=1 in the radiance file header
ncatted -O -h -a near_real_time,global,c,l,1 "$rad_file"

# If the version numbers of the baseline and NRT products may differ,
# then set the version number attribute in the local file copy
current_version_id=$(global_attribute.py --attr version_id $rad_path)
if test $current_version_id -ne $SDPC_NRT_PROCESSING_VERSION ; then
  ncatted -O -h -a version_id,global,m,l,$SDPC_NRT_PROCESSING_VERSION "$rad_file"
fi

# Retrieve the pre-INR metfile from the archive
# and update the LOCALGRANULEID, SHORTNAME, and VERSIONID values.
/bin/cp "$orig_metfile_path" "${rad_basename}.nc.met"
have_version=$(printf "TEMPO_RAD_L1_V%02d" $current_version_id)
want_version=$(printf "TEMPO_RAD_L1_V%02d" $SDPC_NRT_PROCESSING_VERSION)
have_prefix="TEMPO_RAD_L1"
want_prefix="TEMPO_RAD_L1_NRT"
# Note that LOCALGRANULEID gets modified twice here (be careful if you change this!)
sed -i -e "s/${have_version}/${want_version}/g" -e "s/${have_prefix}/${want_prefix}/g" "${rad_basename}.nc.met"
set_metfile_versionid.py --version $SDPC_NRT_PROCESSING_VERSION "${rad_basename}.nc.met"

run_inr_post ${rad_basename}.nc

/bin/cp $irr_file ${irr_basename}.nc
(run_cloud_o4)
create_file_listing  "$cld_o4_basename"

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# To facilitate generation of level 2 products, put a tar notice file
# in $l2_incoming_nrt. Using a notice file instead of the tar
# file itself minimizes data movement and should improve efficiency.

tar_granule_dir_to_dest "$l1_out_dir"

notify_granule_ready "$l2_incoming_nrt"

perform_cleanup
