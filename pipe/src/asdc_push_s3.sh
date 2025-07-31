#! /usr/bin/env bash

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

#set -e
set -u

if test $# -lt 1 ; then
    echo "Usage: $0 Bucket:Bucket_Directory [DBFILE]"
    exit 0
fi

s3_bucket=$1
if test $# -eq 1 ; then
   source_dbfile="$SDPC_ARCHIVE_DBFILE_NRT"
else
   source_dbfile="$2"
fi

# Per-table limit on the number of results from database query
# (limit<=0 means no limit)
: "${SDPC_ASDC_LIMIT:=0}"

if ! test -f "$source_dbfile" ; then
   echo "asdc_push_s3.sh: nonexistent database file: $source_dbfile"
   exit 0
fi

pdr_dbfile="$SDPC_ARCHIVE_DIR/asdc/pdrs_s3.sqlite"
ASDC_TRACK_UPLOADS="asdc_track_uploads.py --dbfile $source_dbfile"

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

do_asdc_s3_upload()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi

  cwd=$(pwd)
  cd $dir

  file_list="files.lis"
  pdr_list="pdrfiles.lis"
  upload_sequence="upload.lis"

  exclude_list="$SDPC_PIPE_DIR/etc/asdc_exclude.csv"

  # make list of new data product files
  $ASDC_TRACK_UPLOADS --limit $SDPC_ASDC_LIMIT --order desc --list new > ${file_list}.new

  # Convert .met files to .cmr.json
  convert_odlmet_to_cmrjson.py --filter ${file_list}.new --output $file_list

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
  $ASDC_TRACK_UPLOADS --set pending $file_list

  # generate manifest files and upload sequence
  asdc_mkscript.sl --bucket $s3_bucket --pdr $pdr_list --output $upload_sequence $file_list

  if test -f "$SDPC_ASDC_TRANSFER_DISABLE" ; then
     echo "asdc_push.sh: transfer disabled ($SDPC_ASDC_TRANSFER_DISABLE exists)"
     return
  fi

  # perform the upload
  error_flag=0
  asdc_s3.py --bucket $s3_bucket --put $upload_sequence || error_flag=1
  if test $error_flag -ne 0 ; then
     echo "*** ERROR: upload failed (see $dir)"
     $ASDC_TRACK_UPLOADS --set new $file_list
     return
  fi

  # mark the file status as "uploaded" and record the upload time
  $ASDC_TRACK_UPLOADS --set uploaded $file_list

  # track PDR files for asdc_pull.sh
  asdc_files.py --dbfile $pdr_dbfile --add $(cat $pdr_list)
}

num=$($ASDC_TRACK_UPLOADS --num new)
if test x"$num" = x0 ; then
   echo "asdc_push_s3.sh: ASDC ingest status: new: $num"
   exit 0
fi

upload_dir_path="${SDPC_ARCHIVE_DIR}/asdc/$(date -u +%Y/%j/push/tempo_pdr_s3_%Y%jT%H%M%SZ)"

do_asdc_s3_upload $upload_dir_path

num_new=$($ASDC_TRACK_UPLOADS --num new)
num_pending=$($ASDC_TRACK_UPLOADS --num pending)
num_uploaded=$($ASDC_TRACK_UPLOADS --num uploaded)
echo "asdc_push_s3.sh: ASDC ingest status: new: $num_new  pending: $num_pending uploaded: $num_uploaded"

