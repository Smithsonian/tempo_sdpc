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
#         rad_filename = basename of the L1 radiance file
#
# 2. Data products produced are automatically archived.
#    On error, a tar file is stored in repro/L2
#
set -e
set -u

: "${SDPC_PIPE_ID:?SDPC_PIPE_ID not set}"

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
   echo "${PROGNAME}[$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
   exit 1
}

trap error_exit ERR

log_message()
{
   printf "${PROGNAME}[$$]: $1\n"
}

test -d "$SDPC_ROOT" || error_exit "$LINENO: cannot access SDPC_ROOT directory: $SDPC_ROOT"

# initialize before loading tar notice file
radref_file=""

# Sourcing the tar file notice defines the variables:
# tar_host = machine with tar file on local disk
# tar_host_file_path = path to tar file on $tar_host
# granule_arch_dir_path = path to L2 archive directory for this granule
# rad_filename = basename of the L1 radiance file
# Optionally: radref_file = basename of radiance reference file
# Optionally: redefine SDPC_LEVEL2_PRODUCTS
. $tar_file_notice

tar_file_basename="$(basename $tar_host_file_path)"
# Trim any basename characters following and including '.',
# but ignoring '.' when it's the first character:
tar_file_basename_sans_ext="${tar_file_basename%.*}"
: "${SDPC_GRANULE_LABEL:=$tar_file_basename_sans_ext}"
export SDPC_GRANULE_LABEL

# ensure upper case product list tokens
level2_products="$SDPC_LEVEL2_PRODUCTS"
level2_products=${level2_products^^}

product_list_tokens="$(echo $level2_products | tr , ' ')"

case "$product_list_tokens" in
   *O3PROF*)
      have_o3p="yes"
      product_list_sans_o3p="$(echo $product_list_tokens | sed -e s/O3PROF//)"
      ;;
   *)
      have_o3p=""
      product_list_sans_o3p="$product_list_tokens"
esac

# Some products may need to wait for a radiance reference file:
if test $SDPC_RADREF_ENABLE -ne 0 && test -n "$SDPC_RADREF_PRODUCTS" ; then

   if test -z "$radref_file" ; then
      # We weren't given a radref filename, but we can try to search for one.
      # If the search fails, radref_file will be the empty string, and we'll
      # try the next alternative.
      if test $SDPC_RADREF_SEARCH -ne 0 ; then
         # When "search" is enabled, then we're allowed a relatively large search window
         radref_file=$(select_radref.py $rad_filename)
      else
         # When "search" is disabled, we can still check to see if a radref for
         # the current scan has already been generated.
         radref_file=$(select_radref.py --thisscan $rad_filename)
      fi
      if test -n "$radref_file" && test -f $radref_file ; then
         printf "radref_file=\"$radref_file\"\n" >> $tar_file_notice
      fi
   fi

   if test -n "$radref_file" ; then
      # radref_file is non-empty.
      # If we've been given the full radref path, and the file is actually there,
      # then we're good to go.  Otherwise, we've presumably been given the basename,
      # and we need to search the archive to get the full path.  If the archive has
      # not yet registered the file, then we may need to wait for it.
      if ! test -f $radref_file ; then
         while true ; do
             # First, check if the RADREF_L1 table exists:
             have_radref=$(sqlite3 -readonly -cmd ".timeout 10000" $SDPC_ARCHIVE_DBFILE "select count(*) from sqlite_master where type='table' and name='RADREF_L1';")
             if test "$have_radref" -ne 0 ; then
                # Ok, the table exists, so we can look for the file:
                radref_path=$(sqlite3 -readonly -cmd ".timeout 10000" $SDPC_ARCHIVE_DBFILE "select path from RADREF_L1 where filename=\"$radref_file\";")
                if test -n "$radref_path" ; then
                   sed -i -e "s,radref_file=$radref_file,radref_file=$radref_path," $tar_file_notice
                   break
                fi
             fi
             sleep 30
         done
      fi
   else
      # We don't have a radref file.
      # If we're generating products that need it, then we will
      # delay generating them until later when the radref file is ready:

      products_that_must_wait=""
      products_that_can_proceed=""
      for p in $product_list_sans_o3p ; do
          case "$SDPC_RADREF_PRODUCTS" in
             *$p*)
                products_that_must_wait="$products_that_must_wait $p"
                ;;
             *)
                products_that_can_proceed="$products_that_can_proceed $p"
                ;;
          esac
      done

      product_list_sans_o3p="$products_that_can_proceed"

      if test -n "$products_that_must_wait" ; then
         # We are generating products that need to wait until a radref becomes available.
         # On $tar_host, create a new hardlink to the tar file to preserve it
         # until later when the radref becomes available:
         tar_host_file_path_wait="${tar_host_file_path}_radref"
         ssh $tar_host ln $tar_host_file_path $tar_host_file_path_wait
         # Create a new tar notice file that refers to the new hardlink we just created.
         radref_wait_dir="$SDPC_PIPE_DIR/stage/granules/level2_input/radref_pending"
         mkdir -p $radref_wait_dir
         basename_sans_extname="$(basename $tar_file_notice .tar | sed -e s/^.//)"
         # change the notice file basename to avoid filename conflicts
         tar_file_notice_wait="$radref_wait_dir/${basename_sans_extname}_radref.tar"
         printf "tar_host=$tar_host\n" > $tar_file_notice_wait
         printf "tar_host_file_path=$tar_host_file_path_wait\n" >> $tar_file_notice_wait
         printf "granule_arch_dir_path=$granule_arch_dir_path\n" >> $tar_file_notice_wait
         printf "rad_filename=$rad_filename\n" >> $tar_file_notice_wait
         printf "export SDPC_LEVEL2_PRODUCTS=\"$products_that_must_wait\"\n" >> $tar_file_notice_wait
      fi
   fi
fi

if test x"$have_o3p" != x ; then
  # load o3prof config parameters
  . $SDPC_PIPE_DIR/etc/o3prof_config.sh
  # generate a product only if this scene is selected
  have_o3p=$(o3p_select.sl --list $o3p_scan_select $SDPC_GRANULE_LABEL)
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

  o3p_host_list=$(seq 0 $((o3p_num_hosts-1)))

  for k in $o3p_host_list ; do
     # To enable parallel cleanup, make a per-host hard link to the tar file
     tar_host_file_path_alias="${tar_host_file_path}_${k}"
     tar_file_notice_alias="${tar_file_notice}_${k}"
     ssh $tar_host ln $tar_host_file_path $tar_host_file_path_alias
     printf "tar_host=$tar_host\n" > $tar_file_notice_alias
     printf "tar_host_file_path=$tar_host_file_path_alias\n" >> $tar_file_notice_alias
     printf "granule_arch_dir_path=$granule_arch_dir_path\n" >> $tar_file_notice_alias
     printf "rad_filename=$rad_filename\n" >> $tar_file_notice_alias
  done
fi

# Now that a private copy (hard link) of the input data is ready for every batch
# job, it's ok to submit all the batch jobs.  Submit the non-o3p batch jobs
# first, to give them a better chance of running "soon" (e.g. when all
# cluster nodes are busy).

slurm_logdir="$SDPC_PIPE_DIR/log/level2/slurm"

if test x"$product_list_sans_o3p" != x ; then
  num_products=$(echo "$product_list_sans_o3p" | wc -w)
  jid=$(sbatch --job-name=L2 --parsable --quiet \
         --comment=$SDPC_GRANULE_LABEL \
         --chdir $l2_run_dir \
         --nodes=1-1 --ntasks=$num_products --ntasks-per-core=1 \
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
  # Singleton dependency requires a job-name unique to this pipeline.
  job_o3p="O3PROF:${SDPC_GRANULE_LABEL}:${SDPC_PIPE_ID}"

  for k in $o3p_host_list ; do
     host_spec="${k}-${o3p_num_hosts}"
     tar_file_notice_alias="${tar_file_notice}_${k}"

     jid=$(sbatch --job-name="$job_o3p" --parsable --quiet \
            --comment=$SDPC_GRANULE_LABEL \
            --partition="$o3p_partition" \
            --nodes=1-1 --ntasks=$o3p_ntasks_per_host --ntasks-per-core=1 \
            --chdir=$l2_run_dir \
            --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.o3prof_batch-%j.out" \
            o3prof_batch.sh "$host_spec" "$tar_file_notice_alias")
     log_message "submitted sbatch $jid: o3prof_batch.sh: $SDPC_GRANULE_LABEL: O3PROF:$k"

  done

  # When all submitted o3p jobs finish, all the o3p blocks will be in the archive.
  # Any node can perform the merge using the previously constructed path,
  jid=$(sbatch --dependency=singleton --job-name="$job_o3p" --parsable \
         --comment=$SDPC_GRANULE_LABEL --quiet \
         --partition="$o3p_partition" \
         --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.o3prof_merge-%j.out" \
         o3prof_merge.sh $granule_arch_dir_path)
  log_message "submitted sbatch $jid: o3prof_merge.sh: $SDPC_GRANULE_LABEL"

fi
