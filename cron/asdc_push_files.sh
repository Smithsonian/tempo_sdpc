#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"

# Max number of files to upload in one batch.
: "${SDPC_ASDC_LIMIT:=500}"

set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $0 USER@HOST <dbfile-path>"
    exit 0
fi
user_at_host=$1
dbfile=$2

# lftp will need an ssh-agent with the relevant keys loaded:
# The file should contain something like:
# export SSH_AUTH_SOCK=/home/temposdpc/.ssh/ssh-agent
# export SSH_AGENT_PID=26882
agent_env_file="$HOME/.ssh/ssh-agent-env"
if ! test -r $agent_env_file ; then
   echo "*** Error: can't find ssh-agent config file: $agent_env_file"
   exit 1
fi
. $agent_env_file

export PATH="${SDPC_ANCILLARY_ROOT}/bin:$PATH"

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
  lftp -f $script
}

dbfile_name=$(basename $dbfile)

num=$(asdc_files.py --dbfile $dbfile --num new)
echo "$dbfile_name push status: new: $num"
if test x"$num" = x0 ; then
   exit 0
fi

push_dir="${SDPC_ANCILLARY_ROOT}/asdc/${dbfile_name}/push"
upload_dir_path="${push_dir}/$(date -u +%Y/${dbfile_name}_pdr_%Y%jT%H%M%SZ)"
do_asdc_upload $upload_dir_path

# log num new/pending:
num=$(asdc_files.py --dbfile $dbfile --num new)
num_after=$(asdc_files.py --dbfile $dbfile --num pending)
echo "$dbfile_name push status: new: $num  pending: $num_after"

