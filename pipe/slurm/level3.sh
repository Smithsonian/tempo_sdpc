#! /bin/sh

# 0. This script is run by cachemon, upon completion of all level 2
#    data products of a particular type from a single scan.
#    The script is triggered by the appearance of a file in
#      $SDPC_ARCHIVE_DIR/registry/scans
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
   printf "*** Error: ${PROGNAME}: $1"
   exit 1
}

# loading $pathlist_file defines these variables:
# product_name = e.g. HCHO_L2
# l3_path = path to target Level 3 data product to be generated
# l2_paths = space-delimited list of level 2 data product files
. $pathlist_file

# Run L2_split on NO2_L2 data products

if test x"$product_name" = x"NO2_L2" ; then
   L2_split -c $SDPC_ROOT/etc/l2_split.cfg $l2_paths || error_exit "L2_split failed"
fi

# Run L2_regrid on all L2 data products

l3_target_dir=$(dirname $l3_path)
l3_basename=$(basename $l3_path)

mkdir -p $l3_target_dir || error_exit "mkdir -p $l3_target_dir failed"

regrid_list="$l3_target_dir/TEMPO_${product_name}.lis"
echo $l3_basename > $regrid_list
for f in $l2_paths ; do
   echo $f >> $regrid_list
done

if test x"$product_name" = x"O3PROF_L2" ; then
   l2_regrid_cfg="${SDPC_ROOT}/etc/l3_o3p.cfg"
else
   l2_regrid_cfg="${SDPC_ROOT}/etc/l3.cfg"
fi

(cd $l3_target_dir && L2_regrid $l2_regrid_cfg) || error_exit "L2_regrid failed"

insert_fixed_metadata.py $l3_path

register_symlink="$SDPC_ARCHIVE_DIR/registry/incoming/$l3_basename"
ln -s $l3_path $register_symlink

/bin/rm -f $pathlist_file
