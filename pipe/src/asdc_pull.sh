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
asdc_user=$(echo $user_at_host | cut -d@ -f1)
asdc_host=$(echo $user_at_host | cut -d@ -f2)

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

# The -E option means "delete source file after successful transfer"
emit_script()
{
   script=$1

cat << EOF > $script
open --user $asdc_user --password DUMMY sftp://$asdc_host
set xfer:log-file lftp_pan.log
set xfer:clobber yes
cd ingest/tempo
mget -E TEMPO*.PAN
exit
EOF
}

do_asdc_download()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi

  cd $dir

  script="lftp.script"

  # make lftp script and run it
  emit_script $script
  lftp -f $script

  # If PANs were retrieved, process them:
  panfiles=$(find . -maxdepth 1 -name "TEMPO*.PAN")
  if test x"$panfiles" != x"" ; then
     asdc_track_uploads.py --pans $panfiles
  fi
}

download_dir_path="${SDPC_ARCHIVE_DIR}/asdc/pull/$(date -u +%Y/%j/tempo_pan_%Y%jT%H%M%SZ)"

do_asdc_download $download_dir_path

