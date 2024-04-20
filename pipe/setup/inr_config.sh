#! /bin/sh

: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

assert_dir_exists()
{
   dir=$1

   if ! test -e $dir ; then
      echo "Directory does not exist: $dir"
      exit 1
   fi
}

inr_run_dir="$SDPC_PIPE_DIR/inr"

INR_CONFIG_SRCDIR="$SDPC_ROOT/etc/inr"
assert_dir_exists $INR_CONFIG_SRCDIR

INR_CONFIG_TARGET_DIR="$inr_run_dir/config"
assert_dir_exists $INR_CONFIG_TARGET_DIR

# INR_REFDATA_DIR contains:
#   ellipsoid altitude map
#   GNA polygon
#   JPL ephemeris
INR_REFDATA_DIR="$SDPC_INRSW_ROOT/config"

# INR_IERS_DIR contains:
#   IERS bulletin A
#INR_IERS_DIR="$SDPC_INRSW_ROOT/config"
INR_IERS_DIR="$SDPC_ANCILLARY_ROOT/var/iers/files"

LOGGING_CONF_FILE="TempoLogging.conf"
PIPELINE_CONF_FILE="TempoPipelineInterfaceSAO.conf"

# INR logging
sed -e s,@INR_PROCESSING_ROOT@,$inr_run_dir,g \
       $INR_CONFIG_SRCDIR/${LOGGING_CONF_FILE}.in > $INR_CONFIG_TARGET_DIR/$LOGGING_CONF_FILE
# INR pipeline
sed -e s,@INR_PROCESSING_ROOT@,$inr_run_dir,g \
    -e s,@INR_IERS_DIR@,$INR_IERS_DIR,g \
    -e s,@INR_REFDATA_DIR@,$INR_REFDATA_DIR,g \
    -e s,@SDPC_REFDATA_DIR@,$SDPC_REFDATA_DIR,g \
    -e s,@SDPC_ANCILLARY_ROOT@,$SDPC_ANCILLARY_ROOT,g \
   $INR_CONFIG_SRCDIR/${PIPELINE_CONF_FILE}.in > $INR_CONFIG_TARGET_DIR/$PIPELINE_CONF_FILE

# Create a backup copy of the original config file to facilitate cold restart
# [obsolete as of INR SW r3.0.0]
#/bin/cp $INR_CONFIG_TARGET_DIR/${PIPELINE_CONF_FILE} \
#        $INR_CONFIG_TARGET_DIR/${PIPELINE_CONF_FILE}.cold_start
