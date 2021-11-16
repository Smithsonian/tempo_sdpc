#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"

set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $0 USER@HOST <dbfile-path>"
    exit 0
fi
user_at_host=$1
dbfile=$2

export PATH="${SDPC_ANCILLARY_ROOT}/src:$PATH"

do_asdc_upload()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi

  cd $dir
  file_list="files.lis"
  script="lftp.script"

  # make list of new data product files
  asdc_files.py --dbfile $dbfile --list new > $file_list

  # mark the new files as "pending"
  asdc_files.py --dbfile $dbfile --set pending $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --dest $user_at_host --output $script $file_list

  # perform the upload
  lftp -f $script > /dev/null 2>&1
}

dbfile_dir=$(dirname $dbfile)
dbfile_name=$(basename $dbfile)

num=$(asdc_files.py --dbfile $dbfile --num new)
echo "$dbfile_name push status: new: $num"
if test x"$num" = x0 ; then
   exit 0
fi

upload_dir_path="${dbfile_dir}/asdc/push/$(date -u +${dbfile_name}_pdr_%Y%jT%H%M%SZ)"
do_asdc_upload $upload_dir_path

# log num new/pending:
num=$(asdc_files.py --dbfile $dbfile --num new)
num_after=$(asdc_files.py --dbfile $dbfile --num pending)
echo "$dbfile_name push status: new: $num  pending: $num_after"

