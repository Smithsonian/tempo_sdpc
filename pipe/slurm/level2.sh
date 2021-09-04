#! /bin/sh

# 0. This script generates level2 data products by running
#    slurm batch jobs in parallel (level2_batch.sh, o3prof_batch.sh).
#
# 1. Command-line arguments provide the basename of the radiance granule,
#    and a tar-file notice.  The tar-file notice defines the following
#    variables:
#         tar_host = name of the machine where the tar file resides
#         tar_host_file_path = path to the tar file on $tar_host
#         granule_arch_dir_path = archive subdirectory for this granule
#
# 2. Data products produced are automatically archived.
#    On error, a tar file is stored in repro/L2
#
set -e
set -u

if test -z "${SDPC_LEVEL2_PRODUCTS}" ; then
   echo "SDPC_LEVEL2_PRODUCTS is not set"
   exit 1
fi

if test $# -ne 2 ; then
  echo "Usage: $0 <tar-file-notice-path> <l2_run_dir>"
  exit 1
fi

tar_file_notice="$1"
l2_run_dir="$2"

PROGNAME="$(basename $0)"
error_exit()
{
   echo "${PROGNAME}[$$]: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

log_message()
{
   printf "${PROGNAME}[$$]: $1\n"
}

test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# Sourcing the tar file notice defines the variables:
# tar_host = machine with tar file on local disk
# tar_host_file_path = path to tar file on $tar_host
# granule_arch_dir_path = path to L2 archive directory for this granule
. $tar_file_notice

tar_file_basename_sans_extname="$(basename $tar_host_file_path .tar)"
: "${SDPC_GRANULE_LABEL:=$tar_file_basename_sans_extname}"
export SDPC_GRANULE_LABEL

# ensure upper case product list tokens
level2_products="$SDPC_LEVEL2_PRODUCTS"
level2_products=${level2_products^^}

product_list_tokens="$(echo $level2_products | tr , ' ')"

have_o3p=""
product_list_sans_o3p=""
for p in $product_list_tokens ; do
   if test x"$p" = x"O3PROF" ; then
      have_o3p="yes"
   else
      product_list_sans_o3p="$p $product_list_sans_o3p"
   fi
done

# When O3PROF is selected, we generate the product only
# when this granule is also a selected scene:
if test x"$have_o3p" != x ; then
: "${SDPC_O3PROF_SCAN_STEP:=3}"
: "${SDPC_O3PROF_SCAN_OFFSET:=0}"
  have_o3p=$(o3p_select.sl --step $SDPC_O3PROF_SCAN_STEP --offset $SDPC_O3PROF_SCAN_OFFSET $SDPC_GRANULE_LABEL)
  if test x"$have_o3p" != x ; then
     log_message "O3PROF selected: $SDPC_GRANULE_LABEL"
  else
     log_message "O3PROF excluded: $SDPC_GRANULE_LABEL"
  fi
fi

# Because the o3p batch jobs go to different compute hosts, we use a hard link
# to provide each o3p batch job with its own private copy of the input data,
# to be deleted upon job completion.
# To avoid a race condition, it's important to set up the input data hard links
# for _all_ of the batch jobs, before _any_ of the batch jobs are submitted.
# This avoids a race condition where, e.g. the set of (hcho,no2,o3t)
# finishes and removes the original tar file before the (o3p) job can
# create hard links to the original tar file and then start running.

if test x"$have_o3p" != x ; then

  # FIXME: put these parameters in a control file somewhere(?)
  ntasks_per_op3_host=20
  num_o3p_hosts=3
  o3p_partition="part2"
  o3p_host_list=$(seq 0 $((num_o3p_hosts-1)))

  for k in $o3p_host_list ; do
     # To enable parallel cleanup, make a per-host hard link to the tar file
     tar_host_file_path_alias="${tar_host_file_path}_${k}"
     tar_file_notice_alias="${tar_file_notice}_${k}"
     ssh $tar_host ln $tar_host_file_path $tar_host_file_path_alias
     printf "tar_host=$tar_host\n" > $tar_file_notice_alias
     printf "tar_host_file_path=$tar_host_file_path_alias\n" >> $tar_file_notice_alias
     printf "granule_arch_dir_path=$granule_arch_dir_path\n" >> $tar_file_notice_alias
  done
fi

# Now that a private copy (hard link) of the input data is ready for every batch
# job, it's ok to submit all the batch jobs.  Submit the non-o3p batch jobs
# first, to give them a better chance of running "soon" (e.g. when all
# cluster nodes are busy).

slurm_logdir="$SDPC_RUN_DIR_MASTER/log/level2/slurm"

if test x"$product_list_sans_o3p" != x ; then
  jid=$(sbatch --job-name=L2 --parsable \
         --comment=$SDPC_GRANULE_LABEL \
         --chdir $l2_run_dir \
         --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.level2_batch-%j.out" \
         level2_batch.sh "$tar_file_notice" "$product_list_sans_o3p")
  log_message "submitted sbatch $jid: level2_batch.sh: $SDPC_GRANULE_LABEL: $product_list_sans_o3p"
else
  # If o3p is the only product, we no longer need the primary tar file.
  # Remove both the tar file notice, and the tar file itself.
  ssh $tar_host /bin/rm -f $tar_host_file_path
  /bin/rm -f $tar_file_notice
fi

if test x"$have_o3p" != x ; then

  # On each host in $o3p_host_list, we'll issue a batch of O3 profile tasks.
  # When each batch is completed, the output blocks are archived.
  # When all of the batch jobs are finished, we run the final merge step to
  # merge the archived blocks into a single O3 profile data product.
  # The merge step is triggered by a slurm 'singleton' dependency on the
  # job name "$job_o3p", defined here.
  job_o3p="O3PROF:${SDPC_GRANULE_LABEL}"

  for k in $o3p_host_list ; do
     host_spec="${k}-${num_o3p_hosts}"
     tar_file_notice_alias="${tar_file_notice}_${k}"

     jid=$(sbatch --job-name="$job_o3p" --parsable \
            --comment=$SDPC_GRANULE_LABEL \
            --partition="$o3p_partition" \
            --nodes=1-1 --ntasks=$ntasks_per_op3_host --ntasks-per-core=1 \
            --chdir=$l2_run_dir \
            --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.o3prof_batch-%j.out" \
            o3prof_batch.sh "$host_spec" "$tar_file_notice_alias")
     log_message "submitted sbatch $jid: o3prof_batch.sh: $SDPC_GRANULE_LABEL: O3PROF:$k"

  done

  # When all submitted o3p jobs finish, all the o3p blocks will be in the archive.
  # Any node can perform the merge using the previously constructed path,
  jid=$(sbatch --dependency=singleton --job-name="$job_o3p" --parsable \
         --comment=$SDPC_GRANULE_LABEL \
         --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.o3prof_merge-%j.out" \
         o3prof_merge.sh $granule_arch_dir_path)
  log_message "submitted sbatch $jid: o3prof_merge.sh: $SDPC_GRANULE_LABEL"

fi
