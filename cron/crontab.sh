#! /bin/sh

set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT not set}"
: "${SDPC_LOCKDIR:?SDPC_LOCKDIR not set}"

cmieast_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmieast.sqlite"
cmiwest_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmiwest.sqlite"
geoscf_sqlite="$SDPC_ANCILLARY_ROOT/var/geoscf/geoscf.sqlite"
ims_sqlite="$SDPC_ANCILLARY_ROOT/var/ims/ims.sqlite"

sqlite_backup_dir="$SDPC_ANCILLARY_ROOT/var/backup"
num_backups=9

lockfile_goes="$SDPC_LOCKDIR/goes_cron.lock"
lockfile_geoscf="$SDPC_LOCKDIR/geoscf_cron.lock"

exit_status()
{
   printf "$2"
   exit $1
}

if ! test -d "$SDPC_ANCILLARY_ROOT" ; then
   exit_status 1 "*** Error: cannot access directory: $SDPC_ANCILLARY_ROOT"
fi

cd $SDPC_ANCILLARY_ROOT

if test -f "etc/crontab.conf" ; then
   . etc/crontab.conf
else
   exit_status 1 "*** Error: cannot access config file: $SDPC_ANCILLARY_ROOT/etc/crontab.conf"
fi

if ! test -d $SDPC_LOCKDIR ; then
   mkdir -p $SDPC_LOCKDIR || exit_status 1 "*** Error: failed creating directory: $SDPC_LOCKDIR"
   chmod 700 $SDPC_LOCKDIR
fi

if ! test -d $sqlite_backup_dir ; then
   mkdir -p $sqlite_backup_dir || exit_status 1 "*** Error: failed creating directory: $sqlite_backup_dir"
fi

if test -f "$aws_config_file" ; then
   export AWS_CONFIG_FILE="$aws_config_file"
fi

rotate_backups()
{
   name="$1"

   indices=$(seq $(($num_backups-1)) -1 1)

   for n in $indices ; do
      this="${name}.${n}"
      next="${name}.$(($n+1))"
      if test -f $this ; then
         /bin/mv $this $next
      fi
   done

   if test -f $name ; then
      /bin/mv $name ${name}.1
   fi
}

sqlite_backup()
{
   dbpath="$1"
   backup_path="$sqlite_backup_dir/$(basename $dbpath)"
   rotate_backups $backup_path
   sqlite3 -cmd ".timeout 20000" $dbpath ".backup $backup_path"
}

explain_error_status()
{
  case "$1" in
    16 ) _msg="Download process unable to acquire lock on lockfile: $lockfile_goes"
       ;;
    17 ) _msg="Download process unable to acquire lock on lockfile: $lockfile_geoscf"
       ;;
    18 ) _msg="Upload process unable to acquire lock on lockfile: $lockfile_goes"
       ;;
    19 ) _msg="Upload process unable to acquire lock on lockfile: $lockfile_geoscf"
       ;;
    *) _msg=""
      ;;
  esac
  echo "$_msg"
}

export PATH="${SDPC_ANCILLARY_ROOT}/bin:${SDPC_ROOT}/bin:${SDPC_OTS_ROOT}/bin:$PATH"

tbeg=$(date +%s)

_task=$1

case $_task in
   IERS )
   test x"$state_iers" = xon || exit 0
   pull_iers.sh $iers_source_url
   ;;

   IMS )
   test x"$state_ims" = xon || exit 0
   pull_ims.sh $ims_url $ims_dir
   ;;

   GOES )
   test x"$state_goes" = xon || exit 0
   flock -E 16 -n $lockfile_goes pull_goes.sh $pda_service_account
   ;;

   GEOSCF )
   test x"$state_geoscf" = xon || exit 0
   SDPC_GEOSCF_VERSION="$geoscf_version" flock -E 17 -n $lockfile_geoscf pull_geoscf.sh $geoscf_source_url
   ;;

   ASDC_GOES )
   test x"$state_asdc_goes" = xon || exit 0
   ( flock -E 18 -n 9
     if test "$?" -eq 0 ; then
        $asdc_pull_method $asdc_dropbox $cmieast_sqlite CMIEAST
        $asdc_pull_method $asdc_dropbox $cmiwest_sqlite CMIWEST
        $asdc_push_method $asdc_dropbox $cmieast_sqlite
        $asdc_push_method $asdc_dropbox $cmiwest_sqlite
     fi
   ) 9> $lockfile_goes
   ;;

   ASDC_GEOSCF )
   test x"$state_asdc_geoscf" = xon || exit 0
   ( flock -E 19 -n 9
     if test "$?" -eq 0 ; then
        $asdc_pull_method $asdc_dropbox $geoscf_sqlite GEOSCF
        $asdc_push_method $asdc_dropbox $geoscf_sqlite
     fi
   ) 9> $lockfile_geoscf
   ;;

   ASDC_IMS )
   test x"$state_asdc_ims" = xon || exit 0
   $asdc_pull_method $asdc_dropbox $ims_sqlite IMS
   $asdc_push_method $asdc_dropbox $ims_sqlite
   ;;

   CLEANUP )
   test -f bin/cleanup.sh || exit 0
   bin/cleanup.sh
   ;;

   BACKUP )
   sqlite_backup $cmieast_sqlite
   sqlite_backup $cmiwest_sqlite
   sqlite_backup $geoscf_sqlite
   sqlite_backup $ims_sqlite
   ;;

   * )
   printf "*** Error: unsupported cron task: $_task"
   false
   ;;
esac

exit_status="$?"
tend=$(date +%s)
tdelta=$((tend-tbeg))
msg="$(date -u +%Y%m%d%H%M%SZ): crontab.sh $_task: exit status ${exit_status}: $tdelta seconds"
explain=$(explain_error_status $exit_status)
if test -n "$explain" ; then
   msg="${msg}\n${explain}"
fi
printf "${msg}\n"

if test x"$exit_status" != x0 ; then
   (export _task ; export msg ; envsubst < $SDPC_ANCILLARY_ROOT/etc/alert_message.tmpl | mailx -t)
fi
