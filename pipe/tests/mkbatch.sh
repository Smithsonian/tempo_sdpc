#! /bin/sh

current_file_name="$1"
num_files_expected="$2"

dir=`dirname "$current_file_name"`
num_files=`ls -1 "$dir" | wc -l`

if test "$num_files" -eq "$num_files_expected" ; then
   /bin/mv "$dir"/* .
   /bin/rmdir "$dir"
fi
