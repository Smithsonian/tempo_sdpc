#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -lt 2 ; then
  echo "Usage: $0 <tar-file> <tar-dir>"
  exit 1
fi

tar_file="$1"
tar_dir="$2"

# When the pipeline completes, the processing directory should be empty.
# Delete it, and then delete the original tar file:
/bin/rm $tar_file
/bin/rmdir $tar_dir

if test $# -eq 3 ; then
   tar_unpack_dir="$3"
   if test -d "$tar_unpack_dir" ; then
     /bin/rmdir "$tar_unpack_dir"
   fi
fi
