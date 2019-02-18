#! /bin/sh
#SBATCH --cpus-per-task=4
#SBATCH --output=/dev/null

# Assume this script is started in a writeable directory
# with the path to a tar file name provided
# on the command line.

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
  echo "Usage: $0 <tar-file> <product-list>"
  exit 1
fi

tar_file="$1"
product_list_arg="$2"

test -r $tar_file || exit 1
tar_file_dir=$(basename $tar_file .tar)

# Setup paths to scripts, config files
#
export PATH="$SDPC_ROOT/bin:$PATH"
etc_dir="$SDPC_ROOT/etc"

# unpack the tar file, then move into the new directory:
tar_unpack_dir=$(pwd)
tar xf $tar_file
cd $tar_file_dir

# 1. At this point, we're in a directory containing:
#          * granule_ident.csv
#          * files.lis, containing the file names:
#            - geolocated level 1 radiance file
#            - level 1 irradiance file
#            - level 2 cloud product
#          * the files named in files.lis
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
. ./files.lis
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
  /bin/cp granule_ident.csv files.lis ${rad_basename}.lis $dir
}

remove_redundant_files()
{
   /bin/rm ${rad_basename}.nc ${irr_basename}.nc $cld_file
   /bin/rm files.lis granule_ident.csv ${rad_basename}.lis
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

do_hcho()
{
  job_hcho="hcho:${SDPC_GRANULE_LABEL}"
  jid_hcho=$(sbatch -w $SLURMD_NODENAME --parsable --job-name=$job_hcho run_trace_gas.sh hcho)
  update_job_list $jid_hcho
}

do_no2()
{
  job_no2="no2:${SDPC_GRANULE_LABEL}"
  jid_no2=$(sbatch -w $SLURMD_NODENAME --parsable --job-name=$job_no2 run_trace_gas.sh no2)
  update_job_list $jid_no2
}

do_o3t()
{
  job_o3t="o3t:${SDPC_GRANULE_LABEL}"
  jid_o3t=$(sbatch -w $SLURMD_NODENAME --parsable --job-name=$job_o3t run_o3t.sh)
  update_job_list $jid_o3t
}

product_list="$(echo $product_list_arg | tr -s , ' ')"
if test x"$product_list" = x ; then
   exit 0
fi

for d in $product_list ; do
  init_product_dir $d
done
remove_redundant_files

for prod in $product_list ; do
  do_${prod}
done

if test X"$jid_list" != X ; then
   job_clean="L2-end:$SDPC_GRANULE_LABEL"
   sbatch -w $SLURMD_NODENAME --job-name=$job_clean \
          --dependency=afterany:$jid_list \
          run_L2_cleanup.sh $tar_file $tar_unpack_dir/$tar_file_dir
fi
