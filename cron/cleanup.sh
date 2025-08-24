#! /usr/bin/env bash

PGMNAME="$(basename $0)"

if test $# -ne 1 ; then
   echo "Usage: $PGMNAME rootdir"
   exit 0
fi

rootdir="$1"

if ! test -d $rootdir ; then
   echo "$PGMNAME: nonexistent directory: $rootdir"
   exit 1
fi

replace_old_subdirs_with_tarfiles()
{
  mmin_arg="$1"
  path="$2"

  dirlist=$(find $path -mindepth 2 -maxdepth 2 -type d -name "???" -mmin $mmin_arg)
  if test -z "$dirlist" ; then
     return
  fi

  echo "$PGMNAME: tar old subdirs: $path"

  for dir in $dirlist ; do
     parent_dir="$(dirname $dir)"
     subdir="$(basename $dir)"
     tar czf "${dir}.tar.gz" --remove-files -C "$parent_dir" "$subdir"
  done
}

_asdc_dirs="cmieast.sqlite cmiwest.sqlite geoscf.sqlite ims.sqlite"

for d in $_asdc_dirs ; do
   replace_old_subdirs_with_tarfiles "+2880" "$rootdir/$d"
done
