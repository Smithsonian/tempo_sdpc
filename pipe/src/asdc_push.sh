#! /usr/bin/env bash

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

set -e
set -u

if test $# -ne 1 ; then
    echo "Usage: $0 USER@HOST:dirpath"
    exit 0
fi

user_at_host=$1

# Per-table limit on the number of results from database query
# (limit<=0 means no limit)
: "${SDPC_ASDC_LIMIT:=0}"

pdr_dbfile="$SDPC_ARCHIVE_DIR/asdc/pdrs.sqlite"

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

  cwd=$(pwd)
  cd $dir

  file_list="files.lis"
  pdr_list="pdrfiles.lis"
  script="lftp.script"

  exclude_list="$SDPC_PIPE_DIR/etc/asdc_exclude.csv"

  # make list of new data product files
  asdc_track_uploads.py --limit $SDPC_ASDC_LIMIT --order desc --list new > $file_list

  # apply the upload filter
  if test -f "$exclude_list" ; then
     asdc_exclude_filter.py --filter $exclude_list $file_list
     # If nothing remains, then clean up and return
     if ! test -s $file_list ; then
        /bin/rm -f $file_list ${file_list}.orig
        cd $cwd
        rmdir $dir
        return
     fi
  fi

  # mark the new files as "pending"
  asdc_track_uploads.py --set pending $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --dest $user_at_host --pdr $pdr_list --output $script $file_list

  # perform the upload
  error_flag=0
  lftp -f $script || error_flag=1
  if test $error_flag -ne 0 ; then
     echo "*** ERROR: upload failed (see $dir)"
     asdc_track_uploads.py --set new $file_list
     return
  fi

  # mark the file status as "uploaded" and record the upload time
  asdc_track_uploads.py --set uploaded $file_list

  # track PDR files for asdc_pull.sh
  asdc_files.py --dbfile $pdr_dbfile --add $(cat $pdr_list)
}

num=$(asdc_track_uploads.py --num new)
if test x"$num" = x0 ; then
   echo "asdc_push.sh: ASDC ingest status: new: $num"
   exit 0
fi

upload_dir_path="${SDPC_ARCHIVE_DIR}/asdc/$(date -u +%Y/%j/push/tempo_pdr_%Y%jT%H%M%SZ)"

do_asdc_upload $upload_dir_path

num_new=$(asdc_track_uploads.py --num new)
num_pending=$(asdc_track_uploads.py --num pending)
num_uploaded=$(asdc_track_uploads.py --num uploaded)
echo "asdc_push.sh: ASDC ingest status: new: $num_new  pending: $num_pending uploaded: $num_uploaded"

