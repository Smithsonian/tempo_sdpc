#! /bin/sh

set -u

: "${SDPC_ASDC_DROPBOX:=jhouck@waps.cfa.harvard.edu}"

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"
: "${SDPC_ROOT:?SDPC_ROOT not set}"
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT not set}"

if ! test -d "$SDPC_ANCILLARY_ROOT" ; then
   printf "*** Error: cannot access directory: $SDPC_ANCILLARY_ROOT"
   exit 1
fi

cd $SDPC_ANCILLARY_ROOT

geoscf_sqlite="$SDPC_ANCILLARY_ROOT/geoscf/geoscf.sqlite"
cmig16_sqlite="$SDPC_ANCILLARY_ROOT/goes/cmig16.sqlite"
cmig17_sqlite="$SDPC_ANCILLARY_ROOT/goes/cmig17.sqlite"

export PATH="${SDPC_ANCILLARY_ROOT}/src:${SDPC_ROOT}/bin:${SDPC_OTS_ROOT}/bin:$PATH"

tstamp=$(date -u +%Y%m%d%H%M%SZ)
tbeg=$(date +%s)

_task=$1

case $_task in
   MET )
   dl_met.sl
   ;;

   SNOW )
   pull_ims.sh
   ;;

   GOES )
   pull_goes.sh
   ;;

   GOES_ASDC )
   asdc_pull.sh $SDPC_ASDC_DROPBOX CMIG16 $cmig16_sqlite
   asdc_pull.sh $SDPC_ASDC_DROPBOX CMIG17 $cmig17_sqlite
   asdc_push.sh $SDPC_ASDC_DROPBOX $cmig16_sqlite
   asdc_push.sh $SDPC_ASDC_DROPBOX $cmig17_sqlite
   ;;

   GEOSCF )
   pull_geoscf.sh
   ;;

   GEOSCF_ASDC )
   asdc_pull.sh $SDPC_ASDC_DROPBOX GEOSCF $geoscf_sqlite
   asdc_push.sh $SDPC_ASDC_DROPBOX $geoscf_sqlite
   ;;

   CLEANUP )
   if test -f src/cleanup.sh ; then
      src/cleanup.sh
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
