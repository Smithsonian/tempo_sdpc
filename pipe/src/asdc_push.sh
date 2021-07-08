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

  exclude_list="$SDPC_RUN_DIR_MASTER/etc/asdc_exclude.lis"

  # make list of new data product files
  asdc_track_uploads.py --list new > $file_list

  # apply the upload filter, if any
  if test -f "$exclude_list" ; then
     /bin/cp "$file_list" "${file_list}.orig"
     grep -v '^#' "$exclude_list" | while read -r line
     do
       sed -i "/$line/d" $file_list
     done
     if ! test -s $file_list ; then
        return
     fi
  fi

  # mark the new files as "pending"
  asdc_track_uploads.py --set pending $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --dest $user_at_host --output $script $file_list

  # perform the upload
  lftp -f $script > /dev/null 2>&1
}

num=$(asdc_track_uploads.py --num new)
if test x"$num" = x0 ; then
   echo "asdc_push.sh: ASDC ingest status: new: $num"
   exit 0
fi

upload_dir_path="${SDPC_ARCHIVE_DIR}/asdc/push/$(date -u +%Y/%j/tempo_pdr_%Y%jT%H%M%SZ)"

do_asdc_upload $upload_dir_path

num_after=$(asdc_track_uploads.py --num pending)
echo "asdc_push.sh: ASDC ingest status: new: $num  pending: $num_after"

