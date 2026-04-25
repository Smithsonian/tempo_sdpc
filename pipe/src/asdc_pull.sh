#! /usr/bin/env bash

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

#set -e
set -u

if test -f "$SDPC_ASDC_TRANSFER_DISABLE" ; then
   echo "asdc_pull.sh: transfer disabled ($SDPC_ASDC_TRANSFER_DISABLE exists)"
   exit 0
fi

if test $# -ne 1 ; then
    echo "Usage: $0 USER@HOST:dirpath"
    exit 0
fi

user_at_host=$1
asdc_user=$(echo $user_at_host | cut -d@ -f1)
asdc_host_dirpath=$(echo $user_at_host | cut -d@ -f2)
asdc_host=$(echo $asdc_host_dirpath | cut -d: -f1)
asdc_dirpath=$(echo $asdc_host_dirpath | cut -d: -f2)

pdr_dbfile="$SDPC_ARCHIVE_DIR/asdc/pdrs.sqlite"

script="lftp.script"
remote_pan_list="pan.lis.remote"
pan_list="pan.lis"

if ! test -f "$SDPC_ARCHIVE_DBFILE" ; then
   echo "asdc_pull.sh: nonexistent database file: $SDPC_ARCHIVE_DBFILE"
   exit 0
fi

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

   # Download a list of PAN files.
   # If none exist, there's nothing else to do.

   lftp <<- EOF > $remote_pan_list
	open --user $asdc_user --password DUMMY sftp://$asdc_host
	cd $asdc_dirpath || exit
	glob echo *.PAN
	exit
	EOF

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
      pdr_file=$(echo $pan_file | sed -e s,.PAN,.PDR,)
      pdr_file_status=$(asdc_files.py --dbfile $pdr_dbfile --status $pdr_file)
      if test x"$pdr_file_status" == x0 ; then
         echo "$pan_file" >> $pan_list
      fi
   done

   if ! test -s $pan_list ; then
      return 1
   fi

   # Generate a script to download the filtered PAN list:

   cat <<- EOF > $script
	open --user $asdc_user --password DUMMY sftp://$asdc_host
	set xfer:log-file lftp_pan.log
	set xfer:clobber yes
	cd $asdc_dirpath || exit
	EOF

  cat $pan_list |
  while read -r pan_file
  do
     echo "get -E $pan_file" >> $script
  done

  echo exit >> $script
}

cleanup()
{
  dir=$1
  if ! test -d $dir ; then
     return
  fi
  file_list="$script lftp_pan.log $remote_pan_list $pan_list"
  for f in $file_list ; do
     path="$dir/$f"
     if test -f $path ; then
        /bin/rm $path
     fi
  done
  cd /tmp
  /bin/rmdir $dir
}

do_asdc_download()
{
  dir=$1
  if ! test -d $dir ; then
     mkdir -p $dir
  fi
  cd $dir

  emit_script $script
  retval=$?
  if test $retval -ne 0 ; then
     cleanup $dir
     return
  fi

  lftp -f $script

  panfiles=$(find . -maxdepth 1 -name "TEMPO*.PAN")
  if test x"$panfiles" != x"" ; then
     # process the downloaded PAN files
     asdc_track_uploads.py --pdrdbfile $pdr_dbfile --pans $panfiles
  else
     cleanup $dir
  fi
}

num=$(asdc_track_uploads.py --num uploaded)
if test x"$num" = x0 ; then
   echo "asdc_pull.sh: ASDC ingest status: uploaded products:$num"
   exit 0
fi

num_pdr=$(asdc_files.py --dbfile $pdr_dbfile --num new)
if test x"$num_pdr" = x0 ; then
   echo "asdc_pull.sh: ASDC ingest status: pending PDRs:$num_pdr"
   exit 0
fi

download_dir_path="${SDPC_ARCHIVE_DIR}/asdc/$(date -u +%Y/%j/pull/tempo_pan_%Y%jT%H%M%SZ)"

do_asdc_download $download_dir_path

num_uploaded=$(asdc_track_uploads.py --num uploaded)
num_problem=$(asdc_track_uploads.py --num problem)
echo "asdc_pull.sh: ASDC ingest status: uploaded:$num_uploaded  problem:$num_problem"
