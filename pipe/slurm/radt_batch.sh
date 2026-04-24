#! /usr/bin/env bash
#SBATCH --output=/dev/null

# Assume this script is started in a writeable directory
# with the path to a level3 RADT file list provided
# on the command line.

# exit on error
#set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

fail_with_message()
{
   echo "*** Error: $1"
   exit 1
}

catch()
{
  if test "$1" != "0" ; then
    echo "Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# check that paths are valid
test -d $SDPC_ROOT || fail_with_message "invalid path: SDPC_ROOT=$SDPC_ROOT"
test -d $SDPC_NODE_DIR || fail_with_message "invalid path: SDPC_NODE_DIR=$SDPC_NODE_DIR"
test -d $SDPC_ARCHIVE_DIR || fail_with_message "invalid path: SDPC_ARCHIVE_DIR=$SDPC_ARCHIVE_DIR"

if test $# -ne 1 ; then
  echo "Usage: $0 <file-list>"
  exit 1
fi

pathlist_file="$1"

test -r $pathlist_file || fail_with_message "cannot read file list: $pathlist_file"

# Reading this file defines these variables:
# product_name= .e.g RADT_L1
# l3_path=unused
# l2_paths= space-delimited list of RADT_L1 files
. $pathlist_file

if test -z "$l2_paths" ; then
  fail_with_message "empty granule list"
fi

parent_dir=$(pwd)
work_dir=$(basename $pathlist_file .nc | sed -e s"/^[.]//")
mkdir $work_dir || fail_with_message "creating work_dir=$work_dir"
cd $work_dir

staging_dir="staging"
granule_list="granule_list.txt"
log_file="log_radt_scan.txt"

gather_copies_to_dir()
{
   to_dir="$1"
   paths="$2"
   for p in $paths ; do
       /bin/cp $p $to_dir
   done
}
scatter_from_dir()
{
   from_dir="$1"
   paths="$2"
   for p in $paths ; do
       p_old="${p}.old"
       if test -f $p ; then
          /bin/mv $p $p_old
       fi
       b=$(basename $p)
       /bin/mv $from_dir/$b $p
       if test -f $p_old ; then
          /bin/rm $p_old
       fi
   done
}

# Copy the input files to a local directory for processing
mkdir $staging_dir || fail_with_message "creating staging directory: $staging_dir"
gather_copies_to_dir "$staging_dir" "$l2_paths"

# Create a list of local granules to process:
find $staging_dir -type f > $granule_list
viirsdnb_dir="$SDPC_ANCILLARY_ROOT/var/viirsdnb/mosaics"

update_archive()
{
   # archive log file and plots with first granule of scan
   first_granule=$(echo $l2_paths | cut -d' ' -f1)
   first_granule_dir=$(dirname $first_granule)
   /bin/mv $log_file $first_granule_dir
   first_granule_name=$(basename $first_granule .nc)
   if test -d Output ; then
      xfrm=$(printf "s,^Output,%s," $first_granule_name)
      tar c --remove-files -f $first_granule_dir/${first_granule_name}.tar Output --transform="$xfrm"
   fi
   # update archived granules
   scatter_from_dir "$staging_dir" "$l2_paths"
   for p in $l2_paths ; do
      # update bounding polygon metadata, MD5 checksum, etc.
      bounding_polygon --replace $p
      bounding_polygon_sync.py --src $p ${p}.met > /dev/null 2>&1
      insert_fixed_metadata.py $p
      md5sum $p > ${p}.md5
      # registry service can now clear 'defer' status for this file
      ln -s $p $SDPC_ARCHIVE_DIR/registry/undefer
   done
   # done with these:
   /bin/rm -f $granule_list
   /bin/rmdir --ignore-fail-on-non-empty $staging_dir
}

tar_to_repro_dir()
{
   cd $parent_dir
   /bin/rm -rf $work_dir/$staging_dir
   tar c --remove-files -f $SDPC_PIPE_DIR/L1/repro/${work_dir}.tar $work_dir
}

srun --output=$log_file \
     --nodes=1 --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
     process_radt_scan.sh $granule_list $viirsdnb_dir
if test "$?" -eq 0 ; then
   update_archive
else
   tar_to_repro_dir
fi

cd $parent_dir
if test -d $work_dir ; then
   /bin/rmdir --ignore-fail-on-non-empty $work_dir
fi

/bin/rm -f $pathlist_file
