#! /usr/bin/env bash
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
test -d $SDPC_NODE_DIR || exit 1
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
etc_dir="$SDPC_PIPE_DIR/etc"

radref_file=""
# tar_file_notice is a short script that defines the variables
# tar_host = machine with the tar file on its local disk
# tar_host_file_path = path to the tar file on $tar_host
# rad_filename = basename of the L1 radiance file
# Optionally: radref_file = path to radiance reference file
. $tar_file_notice

# Work in a unique subdirectory to avoid collisions
# between parallel jobs
unique_subdir_name="${SLURM_JOB_ID}_$$"
mkdir $unique_subdir_name
cd $unique_subdir_name

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

# Move into the directory from the unpacked tar file
tar_unpack_dir=$(pwd)
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
#          * SDPC_ROOT, SDPC_NODE_DIR, SDPC_ARCHIVE_DIR
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

solcal_file_list="${rad_basename}.solcal"
destripe_file_list="${rad_basename}.destripe"

init_product_dir()
{
   # Create target subdirectory, with hard links to (rad, irr, cld)
  dir=$1
  /bin/mkdir -p $dir
  /bin/ln ./${rad_basename}.nc ./${irr_basename}.nc ./$cld_file $dir
  /bin/cp pge_input_basenames.lis ${rad_basename}.lis $dir
}

assign_solcal_cache_file()
{
  molecule=$1
  dir=$molecule

  if test -s $solcal_file_list ; then
     path=$(grep "_IRR${molecule}_" $solcal_file_list || true)
     if test -f "$path" ; then
        printf "$path\n" > "$dir/$solcal_file_list"
     fi
  fi
}

assign_radref_file()
{
   molecule=$1
   dir=$molecule

   if test "$molecule" = HCHO ; then
      echo $radref_file > $dir/radref_file
   fi
}

assign_destripe_file()
{
  molecule=$1
  dir=$molecule

  if test -s $destripe_file_list ; then
     path=$(grep "_DSTR${molecule}_" $destripe_file_list || true)
     if test -f "$path" ; then
        printf "$path\n" > "$dir/$destripe_file_list"
     fi
  fi
}

remove_redundant_files()
{
   /bin/rm ${rad_basename}.nc ${irr_basename}.nc $cld_file
   /bin/rm pge_input_basenames.lis ${rad_basename}.lis
}

product_list="$(echo $product_list_arg | tr -s , ' ')"
# If we have no radiance reference file,
# then remove HCHO from the product list
if ! test -f "$radref_file" ; then
    product_list=$(echo $product_list | sed -e "s,HCHO,,")
fi
# If the product list is empty, we're done
if test x"$product_list" = x ; then
   exit 0
fi

for prod in $product_list ; do
  init_product_dir $prod
done
remove_redundant_files

case "$rad_basename" in
   *_NRT_* )
      slurm_logdir="$SDPC_PIPE_DIR/log/level2_nrt/slurm"
     ;;
   * )
      slurm_logdir="$SDPC_PIPE_DIR/log/level2/slurm"
     ;;
esac

# Run a background process for each product.
# Each product script runs the product code via srun so that slurm
# can track the individual job steps.
for prod in $product_list ; do
  if test $prod = O3TOT ; then
     o3tot_log="$slurm_logdir/${rad_basename}.o3tot-${SLURM_JOB_ID}.out"
     o3tot.sh > $o3tot_log 2>&1 &
  else
     tracegas_log="$slurm_logdir/${rad_basename}.tracegas-${SLURM_JOB_ID}.out"
     assign_solcal_cache_file $prod
     assign_radref_file $prod
     assign_destripe_file $prod
     tracegas.sh $prod > $tracegas_log 2>&1 &
  fi
done

# wait for background jobs to exit
wait

/bin/rm -f "$solcal_file_list" "$destripe_file_list"

level2_finish.sh $tar_file_notice $tar_unpack_dir/$tar_file_dir "$tar_unpack_dir" > /dev/null
