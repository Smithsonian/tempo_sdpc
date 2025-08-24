#! /usr/bin/env bash

# Max number of files to upload in one batch.
: "${SDPC_ASDC_LIMIT:=500}"

set -e
set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

if test $# -ne 3 ; then
    echo "Usage: $(basename $0) rootdir USER@HOST:dirpath <dbfile-path>"
    exit 0
fi
rootdir=$1
user_at_host=$2
dbfile=$3

# If the database file doesn't exist, there's nothing to push
if ! test -f $dbfile ; then
    exit 0
fi

# lftp will need an ssh-agent with the relevant keys loaded:
# The file should contain something like:
# export SSH_AUTH_SOCK=/home/temposdpc/.ssh/ssh-agent
# export SSH_AGENT_PID=26882
agent_env_file="$HOME/.ssh/sdpc-agent-env"
if ! test -r $agent_env_file ; then
   echo "*** Error: can't find ssh-agent config file: $agent_env_file"
   exit 1
fi
. $agent_env_file

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
  asdc_files.py --dbfile $dbfile --limit $SDPC_ASDC_LIMIT --list new > $file_list

  # mark the new files as "pending"
  asdc_files.py --dbfile $dbfile --set pending $file_list

  # generate manifest files and upload script
  asdc_mkscript.sl --dest $user_at_host --output $script $file_list

  # perform the upload
  error_flag=0
  lftp -f $script || error_flag=1
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

dbfile_dir="$rootdir/${dbfile_name}"
upload_dir_path="${dbfile_dir}/$(date -u +%Y/%j/push/${dbfile_name}_pdr_%Y%jT%H%M%SZ)"
do_asdc_upload $upload_dir_path

# log num new/pending/uploaded:
num=$(asdc_files.py --dbfile $dbfile --num new)
num_pending=$(asdc_files.py --dbfile $dbfile --num pending)
num_uploaded=$(asdc_files.py --dbfile $dbfile --num uploaded)
echo "$(date -u $tstamp_fmt): $dbfile_name push status: new:$num  pending:$num_pending  uploaded:$num_uploaded"
