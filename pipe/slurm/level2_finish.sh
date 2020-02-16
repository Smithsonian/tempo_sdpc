#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -lt 2 ; then
  echo "Usage: $0 <tar-file-notice> <tar-dir>"
  exit 1
fi

tar_file_notice="$1"
tar_dir="$2"

# tar_file_notice is a short script that defines the variables
# tar_host = machine with the tar file on its local disk
# tar_host_file_path = path to the tar file on $tar_host
. $tar_file_notice
this_host=$(uname -n | cut -d. -f1)
# Remove the original tar file
if test x"$tar_host" != x"$this_host" ; then
   ssh $tar_host /bin/rm -f $tar_host_file_path
else
   /bin/rm -f $tar_host_file_path
fi

# When the pipeline completes, the processing directory should be empty.
# except for the archive_subdir file.
# Delete the processing directory, and the original tar file notice:
/bin/rm -f $tar_file_notice
/bin/rm -f $tar_dir/archive_subdir
/bin/rmdir $tar_dir

if test $# -eq 3 ; then
   tar_unpack_dir="$3"
   if test -d "$tar_unpack_dir" ; then
     /bin/rmdir "$tar_unpack_dir"
   fi
fi
