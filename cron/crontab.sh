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
   if test x"$state_iers" = xon ; then
      pull_iers.sh $iers_source_url
   fi
   ;;

   IMS )
   if test x"$state_ims" = xon ; then
      pull_ims.sh $ims_url $ims_dir
   fi
   ;;

   GOES )
   if test x"$state_goes" = xon ; then
      pull_goes.sh $pda_service_account
   fi
   ;;

   GEOSCF )
   if test x"$state_geoscf" = xon ; then
      pull_geoscf.sh $geoscf_source_url
   fi
   ;;

   ASDC_GOES )
   if test x"$state_asdc_goes" = xon ; then
      cmig16_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmig16.sqlite"
      cmig17_sqlite="$SDPC_ANCILLARY_ROOT/var/goes/cmig17.sqlite"

      asdc_pull_ack.sh $asdc_dropbox $cmig16_sqlite CMIG16 
      asdc_pull_ack.sh $asdc_dropbox $cmig17_sqlite CMIG17 
      asdc_push_files.sh $asdc_dropbox $cmig16_sqlite
      asdc_push_files.sh $asdc_dropbox $cmig17_sqlite
   fi
   ;;

   ASDC_GEOSCF )
   if test x"$state_asdc_geoscf" = xon ; then
      geoscf_sqlite="$SDPC_ANCILLARY_ROOT/var/geoscf/geoscf.sqlite"

      asdc_pull_ack.sh $asdc_dropbox $geoscf_sqlite GEOSCF 
      asdc_push_files.sh $asdc_dropbox $geoscf_sqlite
   fi
   ;;

   CLEANUP )
   if test -f bin/cleanup.sh ; then
      bin/cleanup.sh
   fi
   ;;

   * )
   printf "** Unsupported cron task: $_task"
   exit 1
   ;;
esac

exit_status="$?"

tend=$(date +%s)
tdelta=$((tend-tbeg))

printf "${tstamp}: crontab.sh $_task: exit status ${exit_status}: $tdelta seconds\n"
