# -*- sh -*-
# This include file is to be included by the log/run script of each service

exec 2>&1

if test -z "$SDPC_PIPE_DIR" ; then
   echo "$SERVICE logging will not be started:  SDPC_PIPE_DIR is not defined"
   exit 1
fi

case "$SERVICE" in
  level1a | level1b | level2 )
     mkdir -p $SDPC_PIPE_DIR/log/$SERVICE/slurm
     ;;
  *)
     ;;
esac

exec s6-log -p -b s15000000 t $SDPC_PIPE_DIR/log/$SERVICE
