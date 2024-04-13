#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"

set -e
set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

if test $# -ne 3 ; then
    echo "Usage: $0 USER@HOST:dirpath <dbfile-path> <pan-prefix>"
    exit 0
fi

user_at_host=$1
dbfile=$2
prefix=$3

if ! test -f $dbfile ; then
   echo "File does not exist: $dbfile"
   exit 1
fi

# lftp will need an ssh-agent with the relevant keys loaded.
# The file should contain something like:
# export SSH_AUTH_SOCK=/home/temposdpc/.ssh/ssh-agent
# export SSH_AGENT_PID=26882
agent_env_file="$HOME/.ssh/sdpc-agent-env"
if ! test -r $agent_env_file ; then
   echo "*** Error: can't find ssh-agent config file: $agent_env_file"
   exit 1
fi
. $agent_env_file

asdc_user=$(echo $user_at_host | cut -d@ -f1)
asdc_host_dirpath=$(echo $user_at_host | cut -d@ -f2)
asdc_host=$(echo $asdc_host_dirpath | cut -d: -f1)
asdc_dirpath=$(echo $asdc_host_dirpath | cut -d: -f2)

export PATH="${SDPC_ANCILLARY_ROOT}/bin:$PATH"

# The -E option means "delete source file after successful transfer"
emit_script()
{
   script=$1

cat << EOF > $script
open --user $asdc_user --password DUMMY sftp://$asdc_host
set xfer:log-file lftp_pan.log
set xfer:clobber yes
cd $asdc_dirpath
glob --exist ${prefix}_*.PAN || exit 0
mget -E ${prefix}_*.PAN
exit
EOF
}

cleanup()
{
  dir=$1
  cd /tmp
  /bin/rm -f $dir/lftp.script $dir/lftp_pan.log
  /bin/rmdir $dir
}

do_asdc_download()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi
  cd $dir

  script="lftp.script"
  emit_script $script

  # When there are no files, lftp should exit with zero status
  lftp -f $script

  # If PANs were retrieved, process them:
  panfiles=$(find . -maxdepth 1 -name "${prefix}_*.PAN")
  if test x"$panfiles" != x"" ; then
     asdc_files.py --dbfile $dbfile --pans $panfiles
  else
     cleanup $dir
  fi
}

dbfile_name=$(basename $dbfile)

num=$(asdc_files.py --dbfile $dbfile --num uploaded)
if test x"$num" = x0 ; then
   echo "$(date -u $tstamp_fmt): $dbfile_name pull status: uploaded:$num"
   exit 0
fi

dbfile_dir="${SDPC_ANCILLARY_ROOT}/var/asdc/${dbfile_name}"
download_dir_path="${dbfile_dir}/$(date -u +%Y/%j/pull/${dbfile_name}_pan_%Y%jT%H%M%SZ)"
do_asdc_download $download_dir_path

num_uploaded=$(asdc_files.py --dbfile $dbfile --num uploaded)
num_accepted=$(asdc_files.py --dbfile $dbfile --num accepted)
num_problem=$(asdc_files.py --dbfile $dbfile --num problem)
echo "$(date -u $tstamp_fmt): $dbfile_name pull status: uploaded:$num_uploaded  accepted:$num_accepted  problem:$num_problem"
