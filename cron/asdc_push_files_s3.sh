#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"

# Max number of files to upload in one batch.
: "${SDPC_ASDC_LIMIT:=500}"

#set -e
set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

if test $# -ne 2 ; then
    echo "Usage: $0 Bucket:Bucket_Directory <dbfile-path>"
    exit 0
fi
s3_bucket=$1
dbfile=$2

# If the database file doesn't exist, there's nothing to push
if ! test -f $dbfile ; then
   echo "asdc_push_files_s3.sh: nonexistent database file: $dbfile"
   exit 0
fi

export PATH="${SDPC_ANCILLARY_ROOT}/bin:$PATH"

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
  upload_sequence="upload.lis"

  # make list of new data product files
  asdc_files.py --dbfile $dbfile --limit $SDPC_ASDC_LIMIT --list new > $file_list

  # mark the new files as "pending"
  asdc_files.py --dbfile $dbfile --set pending $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --bucket $s3_bucket --output $upload_sequence $file_list

  # perform the upload
  error_flag=0
  asdc_s3.py --bucket $s3_bucket --put $upload_sequence || error_flag=1
  if test $error_flag -ne 0 ; then
     echo "*** ERROR: upload failed (see $dir)"
     asdc_files.py --dbfile $dbfile --set new $file_list
     return
  fi

  # mark the files as "uploaded"
  asdc_files.py --dbfile $dbfile --set uploaded $file_list
}

dbfile_name=$(basename $dbfile)

num=$(asdc_files.py --dbfile $dbfile --num new)
echo "$(date -u $tstamp_fmt): $dbfile_name push status: new: $num"
if test x"$num" = x0 ; then
   exit 0
fi

dbfile_dir="${SDPC_ANCILLARY_ROOT}/var/asdc/${dbfile_name}"
upload_dir_path="${dbfile_dir}/$(date -u +%Y/%j/push/${dbfile_name}_pdr_%Y%jT%H%M%SZ)"
do_asdc_upload $upload_dir_path

# log num new/pending/uploaded:
num=$(asdc_files.py --dbfile $dbfile --num new)
num_pending=$(asdc_files.py --dbfile $dbfile --num pending)
num_uploaded=$(asdc_files.py --dbfile $dbfile --num uploaded)
echo "$(date -u $tstamp_fmt): $dbfile_name push status: new:$num  pending:$num_pending  uploaded:$num_uploaded"
