# -*- sh -*-
# This include file is to be included by the log/run script of each service

exec 2>&1

if test -z "$SDPC_PIPE_DIR" ; then
   echo "$SERVICE logging will not be started:  SDPC_PIPE_DIR is not defined"
   exit 1
fi

service_dir="$SDPC_PIPE_DIR/log/$SERVICE"

# Create the directory here with umask permissions.
# When s6-log creates the directory, its unreadable to other users.
mkdir -p $service_dir

case "$SERVICE" in
  level1a | level1b | level2 | level1b_nrt1 | level1b_nrt2 | level2_nrt )
     mkdir -p $service_dir/slurm
     ;;
  *)
     ;;
esac

exec s6-log -p -b s15000000 t $service_dir
