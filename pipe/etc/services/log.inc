# -*- sh -*-
# This include file is to be included by the log/run script of each service

exec 2>1&

if test -z "$SDPC_RUN_DIR_MASTER" ; then
   echo "$SERVICE logging will not be started:  SDPC_RUN_DIR_MASTER is not defined"
   exit 1
fi

exec s6-log -b n20 s1000000 t $SDPC_RUN_DIR_MASTER/log/$SERVICE
