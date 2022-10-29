#! /bin/sh

: "${SDPC_INRSW_ROOT:?SDPC_INRSW_ROOT not set, run this with sdpcrun.sh}"
: "${SDPC_MOCK_SPEEDUP:=1.0}"

set -u
set -e

if test $# -ne 1 ; then
   echo "Usage:  $0 <mock_feed.cfg>"
   exit 0
fi

INR_MOCK_FEED_CFG=$(realpath $1)

# needs a logging config file:
INR_CONFIG_SRCDIR="$SDPC_ROOT/etc/inr"
MOCK_LOGGING_CONF_FILE="TempoMockLogging.conf"
sed -e s,@INR_PROCESSING_ROOT@,$SDPC_PIPE_DIR/inr,g \
       $INR_CONFIG_SRCDIR/${MOCK_LOGGING_CONF_FILE}.in > $MOCK_LOGGING_CONF_FILE

_pgm="${SDPC_INRSW_ROOT}/bin/TempoMockPipeline"
exec $_pgm -s $SDPC_MOCK_SPEEDUP $INR_MOCK_FEED_CFG > feed.$SDPC_PIPE_NAME.log 2>&1 &
