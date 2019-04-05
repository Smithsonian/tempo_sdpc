#! /bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null

# 1. Assume this script is started in a writeable directory
#    with the path to a (hidden) Level 0 granule file provided
#    on the command line.  The L0 granule may contain dark,
#    irradiance or radiance frames.
#
#    The following environment variables are assumed to be set:
#          * SDPC_ROOT, SDPC_RUN_DIR
#
# 2. When execution is successful...
#
#    On error...
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
  echo "Usage: $0 <target-basename> <list-of-files>"
  exit 1
fi
granule_basename="$1"
file_list_file="$2"

# including this file should define these variables:
#    granule_path
#    dark_file_path
. "$file_list_file"

# Setup paths to scripts, config files
# current directory, output directories
#
export PATH="$SDPC_ROOT/bin:$PATH"
etc_dir="$SDPC_ROOT/etc"

l0_incoming_dir="$SDPC_RUN_DIR/L0/incoming"
l0_out_dir="$SDPC_RUN_DIR/L0/out"
l0_repro_dir="$SDPC_RUN_DIR/L0/repro"

l1_out_dir="$SDPC_RUN_DIR/L1/out"
l1_repro_dir="$SDPC_RUN_DIR/L1/repro"

l2_incoming_dir="$SDPC_RUN_DIR/L2/incoming"

inr_input_cache="$SDPC_RUN_DIR/L1/radiance_inr_staging"

# Make a working directory with a local copy of the granule file.
work_dir=$(basename $granule_basename .nc)
/bin/mkdir "$work_dir"
cd $work_dir
/bin/cp "$granule_path" "$granule_basename"
/bin/cp "$file_list_file" "${granule_basename}.lis"
chmod u+w "$granule_basename"

work_dir_tarfile="${work_dir}.tar"

mkgranule_ident -o granule_ident.csv -v $SDPC_PROCESSING_VERSION $granule_basename

run_dir=$(pwd)
parent_dir=$(dirname "$run_dir")

tar_product_to_dest_dir()
{
   dest_dir=$1
   /bin/mkdir -p $dest_dir

   cd $run_dir
   /bin/rm $granule_basename

   cd $parent_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile $work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile
}

finish()
{
   tar_product_to_dest_dir "$l0_repro_dir"
}
trap finish EXIT ERR

run_l0_ccd()
{
   output_file="$1"
   dark_option="$2"

   /bin/cp ${etc_dir}/l0_ccd.cfg .

   srun --ntasks=1 --output=log_l0_ccd.txt \
   L0_ccd -i $l0_incoming_dir/telem \
          -o $output_file $dark_option \
          $granule_basename
}

. $SDPC_ROOT/bin/run_wavecal.sh

run_inr_prep()
{
   target_file="$1"

   /bin/cp ${etc_dir}/l1_inr_prep.cfg .

   srun --ntasks=1 --output=log_inr_prep.txt \
   L1_inr_prep $target_file
}

case "${granule_basename}" in
  *drk* )
  output_file=$(mkgranule_name -L 0 -p drkt -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file ""
  tar_out_dir="$l0_out_dir"
  archive_level="L0"
  ;;

  *irr* )
  output_file=$(mkgranule_name -L 1 -p irr -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file "-d $dark_file_path"
  run_wavecal $output_file "0-4"
  tar_out_dir="$l1_out_dir"
  archive_level="L1"
  ;;

  *rad* )
  output_file=$(mkgranule_name -L 1 -p rad -v $SDPC_PROCESSING_VERSION $granule_basename)
  run_l0_ccd $output_file "-d $dark_file_path"
  # run radiance wavelength calibration post-INR, so as not to delay INR
  run_inr_prep $output_file

  tar_out_dir="$l1_out_dir"
  archive_level="L1"

  # We'll need the metadata file again for post-INR processing,
  # so we'll optimistically put it in the L2 incoming cache
  # (because that's a known path that's already available to all nodes).
  metadata_file="${output_file}.met"
  if test -f "$metadata_file" ; then
     /bin/cp "$metadata_file" $l2_incoming_dir
  fi

  rad_tmpfile=$inr_input_cache/.${output_file}
  /bin/cp $output_file $rad_tmpfile
  /bin/mv $rad_tmpfile $inr_input_cache/$output_file
  ;;

  * )
  printf "*** Unsupported granule filename pattern\n"
  exit 1
  ;;
esac

trap - EXIT

if test x"$tar_out_dir" != x ; then
   tar_product_to_dest_dir "$tar_out_dir"
   tarfile_path="$tar_out_dir/${work_dir_tarfile}"
   archive.sl --delete -a $SDPC_ARCHIVE_DIR -l $archive_level $tarfile_path
fi

# Assume the initial L0 granule was archived when it was produced,
# so it's ok to delete this copy once the archive.sl process has
# succeeded.
/bin/rm "$granule_path" "$file_list_file"
