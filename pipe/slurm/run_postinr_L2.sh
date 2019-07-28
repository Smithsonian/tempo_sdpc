#! /bin/sh

set -e
set -u

# Note that the radiance file name be "hidden" (may begin with a ".").

if test $# -ne 2 ; then
  echo "Usage: $0 <rad-path> <run_dir>"
  exit 1
fi

rad_path="$1"
run_dir="$2"

PROGNAME="$0"
error_exit()
{
   echo "${PROGNAME}: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

test -r $rad_path || error_exit "$LINENO: cannot access granule: $rad_path"
test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

export PATH="$SDPC_ROOT/bin:$PATH"

# Parse the path to the post-INR radiance file
rad_basename=$(basename "$rad_path" .nc| sed -e s"/.Smoothed$//" -e s"/^[.]//")

# construct granule label string for slurm job names
export SDPC_GRANULE_LABEL="${rad_basename}"

# FIXME - in operations, a cron job should handle this DB update, e.g. 1x per day
# For testing, we'll update it here, because the irr file may be newly created
filedb -c $SDPC_ROOT/etc/filedb.cfg tempo:irr --update

# Generate file list file on master node
irr_file=$(filedb -c $SDPC_ROOT/etc/filedb.cfg tempo:irr --find --header "$rad_path")
snow_file=$(filedb -c $SDPC_ROOT/etc/filedb.cfg snow --find --header "$rad_path")
met_file_path=$(filedb -c $SDPC_ROOT/etc/filedb.cfg met --find --header "$rad_path")
granule_dir=$(dirname "$rad_path")
file_list_file="$granule_dir/${rad_basename}.lis"
cat <<EOF > $file_list_file
rad_path=${rad_path}
irr_file=${irr_file}
snow_file=${snow_file}
met_file_path=${met_file_path}
EOF

# Run the post-INR pipeline to prepare for L2 product generation:
job_prep_l2="L2-pre:${SDPC_GRANULE_LABEL}"
sbatch --wait --job-name=$job_prep_l2 --chdir $run_dir \
        --nodes=1-1 --ntasks=8 \
        prep_L2.sh "${rad_basename}.nc" "$file_list_file"

# The --wait above ensures that we don't get to here until
# this tar file has been created -- and there's no point in
# proceeding further if we don't have it.
tarfile_path="$SDPC_RUN_DIR_MASTER/L2/incoming/${rad_basename}.tar"
if ! test -f "$tarfile_path" ; then
  echo "*** Error: prep_L2.sh failed: $rad_basename"
  exit 1
fi

product_list="$SDPC_LEVEL2_PRODUCTS"
product_list_tokens="$(echo $product_list | tr -s , ' ')"
product_list_sans_o3p="$(echo $product_list_tokens | sed -e 's/O3P//g' | tr -s ' ' ,)"

have_o3p=""
for p in $product_list_tokens ; do
   if test x"$p" = x"O3P" ; then
      have_o3p="yes"
   fi
done

# FIXME: for now, only generate o3p products for one granule
case "$rad_basename" in
   *G01* )
   ;;
   * )
   # uncomment this to skip other granules:
   ## have_o3p=""
   ;;
esac

# Because the o3p array jobs go to different compute hosts, we use a hard link
# to provide each o3p array job with its own private copy of the input data,
# to be deleted upon job completion.
# To avoid a race condition, it's important to set up the input data hard links
# for ALL of the batch jobs, before ANY of the batch jobs are actually submitted.
# This avoids a race condition where, e.g. the set of (hcho,no2,o3t)
# finishes and removes the original tar file before the (o3p) job can
# create hard links to the original tar file and then start running.

if test x"$have_o3p" != x ; then

  # FIXME: put these parameters in a control file somewhere.
  ntasks_per_op3_host=14
  num_o3p_hosts=3;
  o3p_host_list=$(seq 0 $((num_o3p_hosts-1)))

  # To perform the o3p merge later, we'll need to know archive path for
  # this granule.  Create that path now, and save it:
  o3p_target_arch_dir_path="$(run_o3p_merge.sh path $tarfile_path)"
  if test x"$o3p_target_arch_dir_path" = x ; then
     echo "*** Error: constructing target archive directory path for $tarfile_path"
     exit 1
  fi

  for k in $o3p_host_list ; do
     # To enable parallel cleanup, make a per-host hard link to the tar file
     # (we can make these links now, only because prep_L2.sh ran with --wait)
     tarfile_path_alias="${tarfile_path}_${k}"
     if ! test -f $tarfile_path_alias ; then
        ln $tarfile_path $tarfile_path_alias
     fi
  done
fi

# Now that a private copy (hard link) of the input data is ready for every batch
# job, it's ok to submit all the batch jobs.  Submit the non-o3p batch jobs
# first, to give them a better chance of running "soon" (e.g. when all
# cluster nodes are busy).

if test x"$product_list_sans_o3p" != x ; then
  job_run_l2="L2:${SDPC_GRANULE_LABEL}"
  sbatch --job-name=$job_run_l2 \
         --chdir $run_dir \
         run_L2.sh "$tarfile_path" "$product_list_sans_o3p"
else
  # If o3p is the only product, we no longer need the primary tar file
  /bin/rm $tarfile_path
fi

if test x"$have_o3p" != x ; then

  # On each host in $o3p_host_list, we'll issue an array of O3 profile tasks
  # with --job-name="$job_o3p".  When each array job is completed, the output
  # blocks are archived.  When all of those job steps are finished, we run the
  # final merge step to merge the blocks into a single O3 profile data product.
  # The merge step is triggered by a slurm 'singleton' dependency on the
  # job name "$job_o3p", defined here.
  job_o3p="o3p:${SDPC_GRANULE_LABEL}"

  for k in $o3p_host_list ; do
     host_spec="${k}-${num_o3p_hosts}"
     tarfile_path_alias="${tarfile_path}_${k}"

     # Here, the --wait ensures that the tar file has been unpacked on each
     # compute host and all associated o3p batch jobs have been submitted
     # BEFORE the singleton dependency cleanup batch job is submitted.
     # Without this wait, there's a race condition, where some compute jobs are
     # submitted after the singleton, causing some blocks to be omitted from
     # the final data product file.
     sbatch --job-name="$job_o3p" \
            --wait --nodes=1-1 --ntasks=$ntasks_per_op3_host \
            --chdir=$run_dir \
            run_L2_o3p.sh "$host_spec" "$tarfile_path_alias"
  done

  # When all submitted o3p jobs finish, all the o3p blocks will be in the archive.
  # Any node can perform the merge using the previously constructed path,
  sbatch --dependency=singleton --job-name="$job_o3p" \
         run_o3p_merge.sh merge $o3p_target_arch_dir_path
fi
