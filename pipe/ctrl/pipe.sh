#! /bin/sh

set -e
set -u

PID_DIR="/var/tmp/$USER/cachemon"

# sdpc_setup.sh must provide definitions for:
# SDPC_ROOT = root directory where pipeline processing software is installed
# SDPC_RUN_DIR = root directory where pipeline processing takes place
# SDPC_ARCHIVE_DIR = root directory for the mission archive
. ./sdpc_setup.sh

test -d "$SDPC_ROOT" || exit 1
test -d "$SDPC_ARCHIVE_DIR" || exit 1

# SDPC_RUN_DIR need not exist on this machine at this point.
# However, it must be defined, and the value will be used
# in the processing directory path on the compute nodes.
: "${SDPC_RUN_DIR:?SDPC_RUN_DIR not set}"

CACHEMON="${SDPC_ROOT}/bin/cachemon.sl"
test -f $CACHEMON || exit 1

#PRODUCT_LIST="hcho,no2,o3t,o3p"
PRODUCT_LIST="hcho,no2,o3t"
#PRODUCT_LIST="o3p"

__init_table_lookup()
{
   # FIXME - During operations, cron jobs will do this.
   #         For testing purposes, do it here.
   ${SDPC_ROOT}/bin/filedb -c $SDPC_ROOT/etc/filedb.cfg met --update
   ${SDPC_ROOT}/bin/filedb -c $SDPC_ROOT/etc/filedb.cfg snow --update
}

start()
{
  __init_table_lookup

  PIPE_SRCDIR=`pwd`

  $CACHEMON -d --rename "${PIPE_SRCDIR}/cachemon_L0_pre_inr.cfg" -- "${SDPC_RUN_DIR}/L0"
  $CACHEMON -d --rename "${PIPE_SRCDIR}/cachemon_postinr_L2.cfg" -- "${SDPC_RUN_DIR}/L2" "$PRODUCT_LIST"
}

stop()
{
  PID_LIST=$(ls $PID_DIR)
  for pid in $PID_LIST; do
    kill -HUP $pid
  done
}

status()
{
  PID_LIST=$(ls $PID_DIR)
  if test X"$PID_LIST" = X ; then
    printf "stopped\n"
  elif test -x $(which pstree) ; then
    for pid in $PID_LIST ; do
      pstree -p $pid
    done
  else
    ps u --ppid $PID_LIST
  fi
}

case "$1" in
  start)
    start
    printf "Pipeline started\n"
    ;;
  stop)
    stop
    printf "Pipeline stopped\n"
    ;;
  status)
    status
    ;;
  *)
    printf "Usage: $0 (start | stop | status)\n"
    exit 0
    ;;
esac
