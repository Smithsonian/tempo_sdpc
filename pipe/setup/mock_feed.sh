#! /bin/sh

: "${SDPC_INRSW_ROOT:?SDPC_INRSW_ROOT not set, run this with sdpcrun.sh}"

set -u
set -e

if test $# -ne 1 ; then
   echo "Usage:  $0 <mock_feed.cfg>"
   exit 0
fi

INR_MOCK_FEED_CFG=$(realpath $1)

_pgm="${SDPC_INRSW_ROOT}/bin/TempoMockPipeline"
exec $_pgm $INR_MOCK_FEED_CFG > feed.log 2>&1 &
