#! /usr/bin/env bash

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

#set -e
set -u

if test -f "$SDPC_ASDC_TRANSFER_DISABLE" ; then
   echo "asdc_pull_s3.sh: transfer disabled ($SDPC_ASDC_TRANSFER_DISABLE exists)"
   exit 0
fi

if test $# -lt 1 ; then
    echo "Usage: $0 Bucket:Bucket_Directory [<arg>]"
    echo "  When optional argument is present, baseline products are transferred instead of NRT"
    exit 0
fi

s3_bucket=$1
if test $# -eq 1 ; then
   source_dbfile="$SDPC_ARCHIVE_DBFILE_NRT"
   pdr_dbfile="$SDPC_ARCHIVE_DIR/asdc/pdrs_nrt.sqlite"
else
   source_dbfile="$SDPC_ARCHIVE_DBFILE"
   pdr_dbfile="$SDPC_ARCHIVE_DIR/asdc/pdrs.sqlite"
fi

if ! test -f "$source_dbfile" ; then
   echo "asdc_pull_s3.sh: nonexistent database file: $source_dbfile"
   exit 0
fi

ASDC_TRACK_UPLOADS="asdc_track_uploads.py --dbfile $source_dbfile"

remote_pan_list="pan_s3.lis.remote"
pan_list="pan_s3.lis"

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

emit_pan_list()
{
   # Download a list of PAN files.
   # If none exist, there's nothing else to do.

   asdc_s3.py --bucket $s3_bucket --list --pattern "*.PAN" > $remote_pan_list
   if ! test -s $remote_pan_list ; then
      return 1
   fi

   # Given a list of PAN files on the remote site, filter the list
   # to contain only those PAN files that correspond to PDR files
   # in $pdr_dbfile (e.g. those previously uploaded by this pipeline
   # instance). This reduces the likelihood of downloading a PAN file
   # that doesn't belong to us.

   touch $pan_list
   sed -i 's,[[:space:]],\n,g' $remote_pan_list

   cat $remote_pan_list |
   while read -r pan_file
   do
      pdr_file=$(echo $pan_file | sed -e "s,.PAN,.PDR,")
      pdr_file_status=$(asdc_files.py --dbfile $pdr_dbfile --status $pdr_file)
      if test x"$pdr_file_status" == x0 ; then
         echo "$pan_file" >> $pan_list
      fi
   done

   if ! test -s $pan_list ; then
      return 1
   fi
}

cleanup()
{
  dir=$1
  if ! test -d $dir ; then
     return
  fi
  file_list="$remote_pan_list $pan_list"
  for f in $file_list ; do
     path="$dir/$f"
     if test -f $path ; then
        /bin/rm $path
     fi
  done
  cd /tmp
  /bin/rmdir $dir
}

do_asdc_s3_download()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi
  cd $dir

  emit_pan_list
  retval=$?
  if test $retval -ne 0 ; then
     cleanup $dir
     return
  fi

  asdc_s3.py --bucket $s3_bucket --get $pan_list

  panfiles=$(find . -maxdepth 1 -name "TEMPO*.PAN")
  if test x"$panfiles" != x"" ; then
     # process the downloaded SHORTPAN files
     $ASDC_TRACK_UPLOADS --pdrdbfile $pdr_dbfile --pans $panfiles
     # delete the processed pans from the S3 bucket
     asdc_s3.py --bucket $s3_bucket --remove $pan_list
  else
     cleanup $dir
  fi
}

# Attempt a download only when we're expecting something.
num_uploaded=$($ASDC_TRACK_UPLOADS --num uploaded)
num_problem=$($ASDC_TRACK_UPLOADS --num problem)
if test -f $pdr_dbfile ; then
   num_pdr=$(asdc_files.py --dbfile $pdr_dbfile --num new)
else
   num_pdr=0
fi
try_download=$(($num_uploaded + $num_problem + $num_pdr))

if test $try_download -eq 0; then
   #echo "asdc_pull_s3.sh: ASDC ingest status: uploaded:$num_uploaded  problem:$num_problem  pending PDRs:$num_pdr"
   exit 0
fi

uniqify="$(mktemp XXXX)"
download_dir_path="${SDPC_ARCHIVE_DIR}/asdc/$(date -u +%Y/%j/pull/tempo_pan_s3_%Y%jT%H%M%SZ_${uniqify})"

do_asdc_s3_download $download_dir_path

num_uploaded=$($ASDC_TRACK_UPLOADS --num uploaded)
num_accepted=$($ASDC_TRACK_UPLOADS --num accepted)
num_problem=$($ASDC_TRACK_UPLOADS --num problem)
echo "asdc_pull_s3.sh: ASDC ingest status: uploaded:$num_uploaded  accepted:$num_accepted  problem:$num_problem"
