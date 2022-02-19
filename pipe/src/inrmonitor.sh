#! /bin/sh

if test -z "$SDPC_PIPE_NAME" ; then
   printf "*** inrmonitor.sh: SDPC_PIPE_NAME not set\n"
   exit 1
fi

top=$(dirname $0)/..
. $top/etc/sdpc_env.sh

inr_run_dir="$SDPC_RUN_DIR_MASTER/inr"

export JAVA_HOME="$SDPC_OTS_ROOT/jdk1.8.0_91"
export PATH="$JAVA_HOME/bin:$PATH"

export TEMPO_INRSW_HOME="$SDPC_INRSW_ROOT"
export TEMPO_INRSW_GUI_HOME="$inr_run_dir"
export TEMPO_INRSW_CONFIG_DIR="$inr_run_dir/config"

options="-Dprism.order=sw"

java $options -jar $SDPC_INRSW_ROOT/gui/GranuleViewerApp.jar \
                   $TEMPO_INRSW_CONFIG_DIR/TempoPipelineInterfaceSAO.conf > /dev/null 2>&1
