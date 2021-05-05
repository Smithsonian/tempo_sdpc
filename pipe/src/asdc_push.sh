#! /usr/bin/env bash

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

set -e
set -u

if test $# -ne 1 ; then
    echo "Usage: $0 USER@HOST"
    exit 0
fi

user_at_host=$1

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

do_asdc_upload()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi

  cd $dir
  file_list="files.lis"
  script="lftp.script"

  # make list of files to upload
  asdc_track_uploads.py --list new > $file_list

  # mark the new files as "pending"
  asdc_track_uploads.py --define $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --dest $user_at_host --output $script $file_list

  # perform the upload
  lftp -f $script > /dev/null 2>&1
}

# When there are no new files, silently do nothing
num=$(asdc_track_uploads.py --num new)
if test x"$num" = x0 ; then
   exit 0
fi

upload_dir_path="${SDPC_ARCHIVE_DIR}/asdc/push/$(date -u +%Y/%j/tempo_pdr_%Y%jT%H%M%SZ)"

do_asdc_upload $upload_dir_path

