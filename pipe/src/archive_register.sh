#! /bin/sh

set -e

if test $# -ne 1; then
   printf "Usage: $0 FILE\n"
   exit 0
fi

target_file="$1"

register.py $target_file
if test -L "$target_file" ; then
   /bin/rm $target_file
fi
