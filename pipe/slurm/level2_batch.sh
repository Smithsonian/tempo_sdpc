#! /bin/sh -vx
#SBATCH --cpus-per-task=4
#SBATCH --output=/dev/null

# Assume this script is started in a writeable directory
# with the path to a tar file notice provided
# on the command line.

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# check that paths are valid
test -d $SDPC_ROOT || exit 1
test -d $SDPC_RUN_DIR || exit 1
test -d $SDPC_ARCHIVE_DIR || exit 1

if test $# -ne 2 ; then
  echo "Usage: $0 <tar-file-notice> <product-list>"
  exit 1
fi

tar_file_notice="$1"
product_list_arg="$2"

test -r $tar_file_notice || exit 1

# Setup paths to scripts, config files
#
export PATH="$SDPC_ROOT/bin:$PATH"
etc_dir="$SDPC_ROOT/etc"

# tar_file_notice is a short script that defines the variables
# tar_host = machine with the tar file on its local disk
# tar_host_file_path = path to the tar file on $tar_host
. $tar_file_notice

# Retrieve the tar file and unpack it:
tar_file_basename=$(basename $tar_host_file_path)
this_host=$(uname -n | cut -d. -f1)
if test x"$tar_host" != x"$this_host" ; then
   scp -o StrictHostKeyChecking=no $tar_host:$tar_host_file_path .
   tar xf $tar_file_basename
   /bin/rm $tar_file_basename
else
   tar xf $tar_host_file_path
fi

# move into the directory from the unpacked tar file
tar_unpack_dir=$(pwd)
tar_file_dir=$(basename $tar_file_basename .tar)
cd $tar_file_dir

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

jid_list=""
update_job_list()
{
  jid=$1
  if test X"$jid_list" = X ; then
     jid_list=$jid
  else
     jid_list="$jid:$jid_list"
  fi
}

product_list="$(echo $product_list_arg | tr -s , ' ')"
if test x"$product_list" = x ; then
   exit 0
fi

for prod in $product_list ; do
  init_product_dir $prod
done
remove_redundant_files

for prod in $product_list ; do
  job_label_args="--job-name=$prod --comment=$SDPC_GRANULE_LABEL"
  if test $prod = O3TOT ; then
     jid=$(sbatch -w $SLURMD_NODENAME --parsable $job_label_args o3tot.sh)
  else
     jid=$(sbatch -w $SLURMD_NODENAME --parsable $job_label_args tracegas.sh $prod)
  fi
  update_job_list $jid
done

if test X"$jid_list" != X ; then
   sbatch -w $SLURMD_NODENAME --job-name="L2:finish" --comment=$SDPC_GRANULE_LABEL \
          --dependency=afterany:$jid_list \
          level2_finish.sh $tar_file_notice $tar_unpack_dir/$tar_file_dir
fi
