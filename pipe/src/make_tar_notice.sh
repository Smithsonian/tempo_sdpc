#! /bin/sh

set -e
set -u

if [ $# -eq 0 ]; then
cat <<EOF
Usage: make_tar_notice.sh <tar-file-absolute-path> [DIR]
   DIR = [Optional] Directory to receive the output in a file
         with the same filename as the input tar file.
         When this argument is missing, the script writes to stdout.
         When this argument is present, the directory must exist (it will not be created).
         The target filename must not exist (it will not be over-written).
         The target filename is created atomicly, via a rename.
EOF
  exit 1
fi

error_exit() {
   printf "*** Error: make_tar_notice.sh: $1\n"
   exit 1
}
trap 'error_exit "unknown error $? line $LINENO"' ERR

case "$1" in
  /* ) tar_host_file_path="$1"
  ;;

  * ) error_exit "needs an absolute path"
  ;;
esac

if ! test -f "$tar_host_file_path" ; then
   error_exit "cannot access file: $tar_host_file_path"
fi

dir=$(basename $tar_host_file_path .tar)
granule_arch_dir_path=$(tar -x -f $tar_host_file_path --to-stdout $dir/archive_subdir)

this_hostname_sans_domain=$(uname -n | cut -d. -f1)

generate_notice() {
cat <<EOF
tar_host="$this_hostname_sans_domain"
tar_host_file_path="$tar_host_file_path"
granule_arch_dir_path="$granule_arch_dir_path"
EOF
}

notice_text=$(generate_notice)

if test $# -eq 1; then
  echo "$notice_text"
  exit 0
fi

outdir="$2"

if ! test -d "$outdir" || ! test -w "$outdir" ; then
   error_exit "cannot write to directory $outdir"
fi

file_basename=$(basename $tar_host_file_path)
temp_file="$outdir/.${file_basename}"
target_file="$outdir/${file_basename}"

for f in "$temp_file" "$target_file" ; do
   if test -f "$f" ; then
      error_exit "file exists: $f"
   fi
done

# Ensure the final output file appears atomically
echo "$notice_text" > $temp_file
/bin/mv $temp_file $target_file
