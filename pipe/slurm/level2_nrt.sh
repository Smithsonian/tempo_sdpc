#! /bin/sh

# 0. This script generates level2 data products by running
#    slurm batch jobs in parallel (level2_batch.sh).
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
product_list_sans_o3p="$(echo $product_list_tokens | sed -e s/O3PROF//)"

case "$product_list_sans_o3p" in
   *)
    # do nothing
    ;;

   *HCHO*)
     # Support HCHO destriping correction
     hcho_destripe_enable=$(config_setting destripe.HCHO.enable)
     destripe_file=""
     if test $hcho_destripe_enable -ne 0 ; then
        hcho_destripe_apply=$(config_setting destripe.HCHO.apply)
        hcho_destripe_search=$(config_setting destripe.HCHO.search)
        if test $hcho_destripe_apply -ne 0 && test $hcho_destripe_search -ne 0 ; then
           # may return empty string if search fails
           destripe_file=$(select_destripe.py $rad_filename)
        fi
     fi
     printf "destripe_file=\"$destripe_file\"\n" >> $tar_file_notice
     ;;
esac

slurm_logdir="$SDPC_PIPE_DIR/log/level2_nrt/slurm"

if test x"$product_list_sans_o3p" != x ; then
  num_products=$(echo "$product_list_sans_o3p" | wc -w)
  jid=$(sbatch --job-name=nL2 --parsable \
         --comment=$SDPC_GRANULE_LABEL \
         --chdir $l2_run_dir \
         --nodes=1-1 --ntasks=$num_products --ntasks-per-core=1 \
         --output "$slurm_logdir/${SDPC_GRANULE_LABEL}.level2_batch-%j.out" \
         level2_batch.sh "$tar_file_notice" "$product_list_sans_o3p")
  log_message "submitted sbatch $jid: level2_batch.sh: $SDPC_GRANULE_LABEL: $product_list_sans_o3p"
fi
