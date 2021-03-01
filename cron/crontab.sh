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

export PATH="${SDPC_ROOT}/bin:${SDPC_OTS_ROOT}/bin:$PATH"

tstamp=$(date -u +%Y%m%d%H%M%SZ)
tbeg=$(date +%s)

_task=$1

case $_task in
   MET )
   ./src/dl_met.sl
   ;;

   SNOW )
   ./src/dl_snow_ims.sl
   ;;

   GOES )
   lftp -f src/pda_lftp.script
   ;;

   CLEANUP )
   ./src/deleter.sl
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
