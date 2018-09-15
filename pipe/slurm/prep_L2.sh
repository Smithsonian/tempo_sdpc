#! /bin/sh
#SBATCH --cpus-per-task=2
#SBATCH --output=/dev/null

# 1. Assume this script is started in a writeable directory
#    with the path to a (hidden) radiance granule provided
#    on the command line.
#
#    The following environment variables are assumed to be set:
#          * SDPC_ROOT, SDPC_RUN_DIR
#
# 2. Scripts lookup_irr.sh and lookup_snow.sh retrieve the appropriate
#    irradiance and snow data products required for processing.
#
# 3. When execution is successful, the following tar files are created:
#          tarfile = $SDPC_RUN_DIR/L1/out/${rad_basename}.tar
#      tarfile_rad = $SDPC_RUN_DIR/L1/out/${rad_basename}.rad.tar
#      tarfile_cld = $SDPC_RUN_DIR/L2/out/${rad_basename}.cld.tar
#
#    On error, the relevant tar files go to the corresponding 'repro'
#    directory.
#
#    The tarfile contents are:
#         tarfile = (rad1, irr1, cld2, etc) for input to L2 pipeline
#     tarfile_rad = (rad1, etc) for archiving
#     tarfile_cld = cloud processing directory, for archiving
#
#---------------------------------------------------------------------

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_RUN_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-path> <rad-target-file>"
  exit 1
fi
rad_path="$1"
rad_file="$2"

# Setup paths to scripts, config files
# current directory, output directories
#
export PATH="$SDPC_ROOT/bin:$PATH"
etc_dir="$SDPC_ROOT/etc"

l1_out_dir="$SDPC_RUN_DIR/L1/out"
l1_repro_dir="$SDPC_RUN_DIR/L1/repro"
l2_incoming="$SDPC_RUN_DIR/L2/incoming"
l2_out_dir="$SDPC_RUN_DIR/L2/out"
l2_repro_dir="$SDPC_RUN_DIR/L2/repro"

# Make a working directory with a local copy of the radiance file.
work_dir=$(basename "$rad_file" .nc)
/bin/mkdir "$work_dir"
cd $work_dir
/bin/cp "$rad_path" "$rad_file"
chmod u+w "$rad_file"

run_dir=$(pwd)
parent_dir=$(dirname "$run_dir")
granule_dir=$(basename "$run_dir")

# Input files:
#
irr_file=$(lookup_irr.sh "$rad_file")
snow_file=$(lookup_snow.sh "$rad_file")
rad_basename=$(basename "$rad_file" .nc)
irr_basename=$(basename "$irr_file" .nc)

# Define template product file name
#
lev1_file_fmt=$(mkgranule_name -p %s "${rad_basename}.nc")
lev1_base_fmt=$(basename "$lev1_file_fmt" .nc)
cld_basename=$(printf "$lev1_base_fmt" cldrr)
cld_dir="cloud"

processing_version=1

get_tiepoint_file()
{
   # If the basename has a leading ".", remove it
   rad_path_basename_sans_ext=$(basename $rad_path .nc | sed -e s"/^[.]//")
   rad_path_dir=$(dirname $rad_path)

   tiepoint_path="$rad_path_dir/${rad_path_basename_sans_ext}.Internal.nc"

   if test -f "$tiepoint_path" ; then
      tiepoint_file=$(printf "$lev1_file_fmt" tie)
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
trap finish EXIT

tar_l1_radiance_to_dest()
{
   dest_dir=$1
   mkdir -p "$dest_dir"
   cd $parent_dir
   tarfile_rad="${rad_basename}.rad.tar"

   if test x"$tiepoint_file" != x ; then
      tiepoint_file="$granule_dir/$tiepoint_file"
   fi

   tar cf $dest_dir/.${tarfile_rad} \
       $granule_dir/${rad_basename}.nc \
       $granule_dir/granule_ident.csv \
       $granule_dir/log_inr_post.txt $tiepoint_file \
       $granule_dir/log_polcorr.txt
   /bin/mv $dest_dir/.${tarfile_rad} $dest_dir/${tarfile_rad}

   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L1 $dest_dir/${tarfile_rad}

   /bin/rm -f $granule_dir/log_inr_post.txt $tiepoint_file \
              $granule_dir/log_polcorr.txt
}

. $SDPC_ROOT/bin/run_wavecal.sh

run_inr_post()
{
   radiance_file=$1

   # INR post
   srun --ntasks=1 --output=log_inr_post.txt \
    L1_inr_post -c ${etc_dir}/l1_inr_post.cfg \
                -s $snow_file $radiance_file

   run_wavecal $radiance_file "0-9"

   # polarization correction
   srun --ntasks=1 --output=log_polcorr.txt \
    L1_polcorr -c ${etc_dir}/l1_inr_post.cfg $radiance_file

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
  . ${etc_dir}/cloud/cloud.rc

  mkdir -p $cld_dir
  cd $cld_dir

  radiance_file="../${rad_basename}.he4"
  irradiance_file="../${irr_basename}.he4"
  product_file="${cld_basename}.he5"
  template_pcf="${etc_dir}/cloud/default.pcf.in"

  # Edit the PCF file template:
  sed \
      -e s,@refdata_dir@,$refdata_dir,g \
      -e s,@spectra_dir@,$spectra_dir,g \
      -e s,@product_dir@,$product_dir,g \
      -e s,@radiance_file@,$radiance_file,g \
      -e s,@irradiance_file@,$irradiance_file,g \
      -e s,@product_file@,$product_file,g \
      $template_pcf > $pcf_file

  export PGS_PC_INFO_FILE="$pcf_file"
  export PGSMSG="${SDPC_ROOT}/msgs"

  srun --ntasks=1 --output=log_cloud.txt \
    L1_cloud $cmdline_args

  tar_l2_cloud_to_dest "$l2_out_dir"

  cd $granule_dir
  /bin/mv ${cld_dir}/${cld_basename}.nc .
  /bin/rm ${cld_dir}/*
  /bin/rmdir ${cld_dir}
}

create_file_listing()
{
cat << EOF > files.lis
RAD=${rad_basename}.nc
IRR=${irr_basename}.nc
CLD=${cld_basename}.nc
EOF
}

perform_cleanup()
{
   # Delete the original radiance file
   /bin/rm "$rad_path"

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

get_tiepoint_file

mkgranule_ident -o granule_ident.csv -v $processing_version ${rad_basename}.nc
run_inr_post ${rad_basename}.nc

if ! test -f "$irr_file" ; then
  echo "ERROR:  irradiance file not found:  $irr_file"
  exit 1
fi

/bin/cp $irr_file ${irr_basename}.nc
(run_cloud)

create_file_listing

trap - EXIT
tar_granule_dir_to_dest "$l2_incoming"

perform_cleanup
