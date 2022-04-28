#! /bin/sh

set -u

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT not set}"

if ! test -d "$SDPC_ANCILLARY_ROOT" ; then
   printf "*** Error: cannot access directory: $SDPC_ANCILLARY_ROOT"
   exit 1
fi

cd $SDPC_ANCILLARY_ROOT

if test -f "etc/crontab.conf" ; then
   . etc/crontab.conf
else
   printf "*** Error: cannot access config file: $SDPC_ANCILLARY_ROOT/etc/crontab.conf"
   exit 1
fi

export PATH="${SDPC_ANCILLARY_ROOT}/bin:${SDPC_ROOT}/bin:${SDPC_OTS_ROOT}/bin:$PATH"

tstamp=$(date -u +%Y%m%d%H%M%SZ)
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
   pull_goes.sh $pda_service_account
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
printf "${tstamp}: crontab.sh $_task: exit status ${exit_status}: $tdelta seconds\n"
