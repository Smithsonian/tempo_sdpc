#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"

#set -e
set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

if test $# -ne 3 ; then
    echo "Usage: $0 Bucket:Bucket_Directory <dbfile-path> <pan-prefix>"
    exit 0
fi

s3_bucket=$1
dbfile=$2
prefix=$3

if ! test -f $dbfile ; then
   echo "asdc_pull_ack_s3.sh: nonexistent database file: $dbfile"
   exit 1
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

pan_list="pan_s3.lis"

cleanup()
{
  dir=$1
  cd /tmp
  /bin/rm -f $dir/$pan_list
  /bin/rmdir $dir
}

do_asdc_download()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi
  cd $dir

  # Download a list of PAN files with the specified prefix.
  # If none exist, there's nothing else to do.

  asdc_s3.py --bucket $s3_bucket --list --pattern "${prefix}_*.PAN" > $pan_list
  if ! test -s $pan_list ; then
     return 0
  fi
  asdc_s3.py --bucket $s3_bucket --get $pan_list

  panfiles=$(find . -maxdepth 1 -name "${prefix}_*.PAN")
  # If no PANs were retrieved, delete the working directory and return
  if test -z "$panfiles" ; then
     cleanup $dir
     return 0
  fi

  # If PANs were retrieved, process them.
  # If processing is successful, delete the PANs from the dropbox
  asdc_files.py --dbfile $dbfile --pans $panfiles
  if test $? -eq 0 ; then
     asdc_s3.py --bucket $s3_bucket --remove $pan_list
  fi
}

dbfile_name=$(basename $dbfile)

num_uploaded=$(asdc_files.py --dbfile $dbfile --num uploaded)
num_problem=$(asdc_files.py --dbfile $dbfile --num problem)
try_download=$(($num_uploaded + $num_problem))
if test $try_download -eq 0 ; then
   echo "$(date -u $tstamp_fmt): $dbfile_name pull status: uploaded:$num_uploaded problem:$num_problem"
   exit 0
fi

dbfile_dir="${SDPC_ANCILLARY_ROOT}/var/asdc/${dbfile_name}"
download_dir_path="${dbfile_dir}/$(date -u +%Y/%j/pull/${dbfile_name}_pan_%Y%jT%H%M%SZ)"
do_asdc_download $download_dir_path

num_uploaded=$(asdc_files.py --dbfile $dbfile --num uploaded)
num_accepted=$(asdc_files.py --dbfile $dbfile --num accepted)
num_problem=$(asdc_files.py --dbfile $dbfile --num problem)
echo "$(date -u $tstamp_fmt): $dbfile_name pull status: uploaded:$num_uploaded  accepted:$num_accepted  problem:$num_problem"
