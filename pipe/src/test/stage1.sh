#! /bin/sh

current_file_name="$1"
num_files_expected="$2"

in_dir=`dirname "$current_file_name"`
num_files=`ls -1 "$in_dir" | wc -l`

if test "$num_files" -lt "$num_files_expected" ; then
  exit 0
fi

/bin/mv "$in_dir"/* .

this_dir=`pwd`
above_dir=`dirname "$this_dir"`
this_subdir=`basename "$this_dir"`

temp_tarfile="${above_dir}/${this_subdir}.tar_in"
dest_tarfile="${above_dir}/${this_subdir}.tar"
tar cf "$temp_tarfile" -C "$above_dir" "$this_subdir" || exit 1
/bin/mv "$temp_tarfile" "$dest_tarfile" || exit 1

if test -f "$dest_tarfile" ; then
  /bin/rm xxx*.dat
  /bin/rmdir "$this_dir"
  /bin/rmdir "$in_dir"
fi
