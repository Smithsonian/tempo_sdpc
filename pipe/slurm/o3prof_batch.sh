#! /bin/sh
#SBATCH --output=/dev/null

# Assume this script is started in a writeable directory
# with the path to a tar file notice provided
# on the command line.

# exit on error
#set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_RUN_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <host-spec> <tar-file-notice>"
  exit 1
fi

# Setup paths to scripts, config files
#
etc_dir="$SDPC_RUN_DIR_MASTER/etc"

# host_spec is a string of the form k:N indicating that
# this is the kth host from a set of N
host_spec="$1"
this_host=$(echo $host_spec | cut -d- -f 1)

# tar_file_notice_alias will have the form $DIR/basename.tar_${k}
# where k indicates which host will be processing the file.
tar_file_notice_alias="$2"

# sourcing the tar file notice defines the variable tar_host and tar_host_file_path
. $tar_file_notice_alias
tarfile_dir=$(basename $tar_host_file_path ".tar_${this_host}")

# Unpack the tar file, then move into the new directory.
# Do this in a uniquely named subdirectory to avoid collisions
# between parallel jobs
unique_subdir_name="${SLURM_JOB_ID}_$$"
mkdir $unique_subdir_name
cd $unique_subdir_name
tar_unpack_dir=$(pwd)

# Retrieve the tar file and unpack it.
tar_file_basename=$(basename $tar_host_file_path)
this_host=$(uname -n | cut -d. -f1)
if test x"$tar_host" != x"$this_host" ; then
   scp -o StrictHostKeyChecking=no $tar_host:$tar_host_file_path .
   tar xf $tar_file_basename
   /bin/rm $tar_file_basename
else
   tar xf $tar_host_file_path
fi

cd $tarfile_dir

# 1. At this point, we're in a directory containing:
#          * archive_subdir
#          * pge_input_basenames.lis, containing the file names:
#            - geolocated level 1 radiance file
#            - level 1 irradiance file
#            - level 2 cloud product
#          * the files named in pge_input_basenames.lis
#
#    The following environment variables are assumed to be set:
#          * SDPC_ROOT, SDPC_RUN_DIR, SDPC_ARCHIVE_DIR
#
# 2. When execution is completed, one or more tar files will
#    have been moved to the appropriate destination directories,
#    and the working directory will be empty.
#
# 3. It is assumed that the calling process will delete the
#    working directory.

run_dir=$(pwd)

# get input file names
. ./pge_input_basenames.lis
rad_file=$RAD
irr_file=$IRR
cld_file=$CLD

rad_basename=$(basename $rad_file .nc)
irr_basename=$(basename $irr_file .nc)

init_product_dir()
{
   # Create target subdirectory, with hard links to (rad, irr, cld)
  dir=$1
  /bin/mkdir -p $dir
  /bin/ln ./${rad_basename}.nc ./${irr_basename}.nc ./$cld_file $dir
  /bin/cp pge_input_basenames.lis ${rad_basename}.lis $dir
}

remove_redundant_files()
{
   /bin/rm ${rad_basename}.nc ${irr_basename}.nc $cld_file
   /bin/rm pge_input_basenames.lis ${rad_basename}.lis
}

do_O3PROF()
{
   array_bounds=$(o3prof_util.sh init "$host_spec")
   joblog_file="O3PROF/joblog.${SLURM_JOB_ID}.out"

   parallel --delay .2 -j $SLURM_NTASKS --joblog $joblog_file \
            o3prof_block.sh {1} $run_dir ::: $(seq $array_bounds)

   o3prof_util.sh cleanup "$host_spec"
}

init_product_dir O3PROF
remove_redundant_files

do_O3PROF

level2_finish.sh ${tar_file_notice_alias} $tar_unpack_dir/$tarfile_dir "$tar_unpack_dir" > /dev/null
