#! /bin/sh

: "${SDPC_INRSW_ROOT:?SDPC_INRSW_ROOT not set, run this with sdpcrun.sh}"
: "${SDPC_INR_RUN_DIR:?SDPC_INR_RUN_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

if test $# -ne 1 ; then
   echo "Usage:  $0 <mock_feed.cfg>"
   exit 0
fi

INR_MOCK_FEED_CFG=$1

INR_FEED="${SDPC_INRSW_ROOT}/bin/TempoMockPipeline"

/bin/cp $INR_MOCK_FEED_CFG $SDPC_INR_RUN_DIR
cd $SDPC_INR_RUN_DIR
exec $INR_FEED $INR_MOCK_FEED_CFG > feed.log 2>&1 &
