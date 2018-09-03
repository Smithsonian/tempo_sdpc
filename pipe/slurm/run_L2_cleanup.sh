#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -ne 2 ; then
  echo "Usage: $0 <tar-file> <tar-dir>"
  exit 1
fi

tar_file="$1"
tar_dir="$2"

# When the pipeline completes, the processing directory should be empty.
# Delete it, and then delete the original tar file:
/bin/rm $tar_file
/bin/rmdir $tar_dir
