#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

PGMNAME="$(basename $0)"

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

_asdc_root="$SDPC_ANCILLARY_ROOT/var/asdc"
_asdc_dirs="cmieast.sqlite cmiwest.sqlite geoscf.sqlite ims.sqlite"

for d in $_asdc_dirs ; do
   replace_old_subdirs_with_tarfiles "+2880" "$_asdc_root/$d"
done
