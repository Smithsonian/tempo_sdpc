#! /bin/sh

# 0. This script is run by cachemon, upon completion of all level 2
#    data products of a particular type from a single scan.
#    The script is triggered by the appearance of a file in
#      $SDPC_PIPE_DIR/stage/scans
#    The argument passed to the script is the path to this file.

#set -e
set -u

if test -z SDPC_ROOT ; then
   printf "SDPC environment is not set\n"
   exit 1
fi

if test $# -ne 1 ; then
   echo "Usage: $0 <L2-pathlist-file>"
   exit 1
fi

pathlist_file="$1"

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

log_message()
{
   printf "${PROGNAME}: $1\n"
}

error_exit()
{
   printf "*** Error: ${PROGNAME}: $1\n"
   exit 1
}

public_mirror_symlink()
{
   src_paths=$1

   mirror_dir="$SDPC_PIPE_DIR/public_mirror"

   # Do nothing when the mirror directory is absent
   if ! test -d $mirror_dir ; then
      return 0
   fi

   for src in $src_paths ; do
       bn=$(basename $src)
       case $bn in
         *_L2_*)
            level_dir=L2
         ;;
         *_L3_*)
            level_dir=L3
         ;;
       esac
       day=$(level1_info --localday $src)
       target_dir="$mirror_dir/$day/${level_dir}"
       if ! test -d $target_dir ; then
          mkdir -p $target_dir
       fi
       ln -s $src $target_dir/$bn || error_exit "public_mirror_symlink failed: $bn"
   done
}

# Run L2_split on NO2_L2 data products
no2_l2_split()
{
   l3_path="$1"
   l2_paths="$2"
   # put log file in the the first granule directory
   first_granule=$(echo $l2_paths | cut -d' ' -f1)
   logdir=$(dirname $first_granule)
   first_granule_bn="$(basename $first_granule .nc)"
   log_message "strat/trop separation: $first_granule_bn scan"
   L2_split -v -c $SDPC_PIPE_DIR/etc/l2_split.cfg $l2_paths > $logdir/log_split.txt 2>&1 || error_exit "L2_split failed"
   public_mirror_symlink "$l2_paths"
}

change_asdc_status_defer_to_new()
{
   _is_nrt="$1"
   l2_paths="$2"

   if test $_is_nrt -ne 0 ; then
      dbfile="$SDPC_ARCHIVE_DBFILE_NRT"
   else
      dbfile="$SDPC_ARCHIVE_DBFILE"
   fi

   tmpfile=$(mktemp)
   printf "%s\n" $l2_paths > $tmpfile
   asdc_track_uploads.py --dbfile $dbfile --undefer $tmpfile
   if test "$?" -ne 0 ; then
      error_exit "asdc_track_uploads failed: changing $product_name asdc_status defer to new in $dbfile"
   fi
   /bin/rm -f $tmpfile
}

# Import function to generate radiance reference file
. $SDPC_ROOT/bin/make_radref.sh

# Import function to generate destriping correction files
. $SDPC_ROOT/bin/make_destripe.sh

# Import function to generate background correction files
. $SDPC_ROOT/bin/make_background.sh

# loading $pathlist_file defines these variables:
# product_name = e.g. HCHO_L2
# l3_path = path to target Level 3 data product to be generated
# l2_paths = space-delimited list of level 2 data product files
. $pathlist_file

l3_target_dir=$(dirname $l3_path)
l3_basename=$(basename $l3_path)

is_nrt=0
case "$l3_basename" in
    *_NRT_* )
       is_nrt=1
       ;;
    * )
       ;;
esac

# Optional: perform destriping
_destripe_products=$(echo $SDPC_DESTRIPE_TG | tr , ' ')
for p in $_destripe_products ; do
    if test $p = $product_name ; then
       # This is a no-op if destriping has already been done
       destripe_scan "$l2_paths"
    fi
done

# Optional: perform background correction
_bkgcorr_products=$(echo $SDPC_BKGCORR_TG | tr , ' ')
for p in $_bkgcorr_products ; do
    if test $p = $product_name ; then
       # This is a no-op if background correction has already been done
       bkgcorr_scan "$l2_paths"
    fi
done

case "$product_name" in
  CLDO4_L2 )
     if test $SDPC_RADREF_ENABLE -ne 0 && test $is_nrt -eq 0 ; then
        make_radref "$l2_paths"
     fi
     ;;

  NO2_L2 )
     no2_l2_split "$l3_path" "$l2_paths"
     ;;

  * )
     ;;
esac

# At this point, any L2 products with asdc_status='defer'
# should be ready to upload to ASDC, so we set them to asdc_status='new'
change_asdc_status_defer_to_new "$is_nrt" "$l2_paths"

# Run L2_regrid on all L2 data products

mkdir -p $l3_target_dir || error_exit "mkdir -p $l3_target_dir failed"

regrid_list="$l3_target_dir/TEMPO_${product_name}.lis"
echo $l3_basename > $regrid_list
for f in $l2_paths ; do
   echo $f >> $regrid_list
done

case "$product_name" in
  O3PROF_L2 )
      l2_regrid_cfg="${SDPC_PIPE_DIR}/etc/l3_o3p.cfg"
      ;;
  O3TOT_L2 )
      l2_regrid_cfg="${SDPC_PIPE_DIR}/etc/l3_o3t.cfg"
      ;;
  * )
      l2_regrid_cfg="${SDPC_PIPE_DIR}/etc/l3.cfg"
      ;;
esac

log_message "generating L3 product: $l3_basename"
(cd $l3_target_dir && L2_regrid -v $l2_regrid_cfg > log_regrid_${product_name}.txt 2>&1 ) || error_exit "L2_regrid failed"

insert_fixed_metadata.py $l3_path
fix_met_format.py ${l3_path}.met
fix_nrt_shortname.py $l3_path

register_symlink="$SDPC_ARCHIVE_DIR/registry/incoming/$l3_basename"
ln -s $l3_path $register_symlink

public_mirror_symlink $l3_path

/bin/rm -f $pathlist_file
