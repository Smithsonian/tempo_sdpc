# -*- sh -*-
# This include file is to be included by the service finish script.
#
# This finish script is executed every time the run script dies.
# It is passed two arguments:
#  * the exit code from the run script (256 if killed by a signal)
#  * the number of the signal that killed the run script (if killed)
#
# It the finish script exits with 125, then the service will not be restarted.
#
exitval="$1"
sigval="$2"

if [ -z "$SDPC_ROOT" ]
then
  echo "$SERVICE will not be restarted: SDPC_ROOT is not defined" 1>&2
  exit 125
fi

if [ X$exitval != X256 ]
then
   exit 0
fi

echo "`date --iso-8601=seconds` $SERVICE killed by signal $sigval" 1>&2
