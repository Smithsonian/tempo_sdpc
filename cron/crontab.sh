#! /bin/sh

set -u

tstamp_fmt="+%Y%m%d%H%M%SZ"

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT not set}"
: "${SDPC_LOCKDIR:?SDPC_LOCKDIR not set}"

exit_status()
{
   status="$1"
   printf "$2"
   exit $status
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
fi

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
   flock -n $SDPC_LOCKDIR/goes_cron.lock pull_goes.sh $pda_service_account
   ;;

   GEOSCF )
   test x"$state_geoscf" = xon || exit 0
   pull_geoscf.sh $geoscf_source_url
   ;;

   ASDC_GOES )
   test x"$state_asdc_goes" = xon || exit 0
   cmieast_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmieast.sqlite"
   cmiwest_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmiwest.sqlite"
   asdc_pull_ack.sh $asdc_dropbox $cmieast_sqlite CMIEAST
   asdc_pull_ack.sh $asdc_dropbox $cmiwest_sqlite CMIWEST
   asdc_push_files.sh $asdc_dropbox $cmieast_sqlite
   asdc_push_files.sh $asdc_dropbox $cmiwest_sqlite
   ;;

   ASDC_GEOSCF )
   test x"$state_asdc_geoscf" = xon || exit 0
   geoscf_sqlite="$SDPC_ANCILLARY_ROOT/var/geoscf/geoscf.sqlite"
   asdc_pull_ack.sh $asdc_dropbox $geoscf_sqlite GEOSCF
   asdc_push_files.sh $asdc_dropbox $geoscf_sqlite
   ;;

   ASDC_IMS )
   test x"$state_asdc_ims" = xon || exit 0
   ims_sqlite="$SDPC_ANCILLARY_ROOT/var/ims/ims.sqlite"
   asdc_pull_ack.sh $asdc_dropbox $ims_sqlite IMS
   asdc_push_files.sh $asdc_dropbox $ims_sqlite
   ;;

   CLEANUP )
   test -f bin/cleanup.sh || exit 0
   bin/cleanup.sh
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
printf "${msg}\n"

if test x"$exit_status" != x0 ; then
   (export _task ; export msg ; envsubst < $SDPC_ANCILLARY_ROOT/etc/alert_message.tmpl | mailx -t)
fi
