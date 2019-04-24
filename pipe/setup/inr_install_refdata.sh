#! /bin/sh

: "${SDPC_INR_RUN_DIR:?SDPC_INR_RUN_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

# INR_GOES_REFDATA_DIR contains:
#   GOES imagery
#   GOES rsr files
INR_GOES_REFDATA_DIR=/data/tempo/sdpc/inr/synth_granules_v2/FullGranules

ln -s -t $SDPC_INR_RUN_DIR/Staging/rsr $INR_GOES_REFDATA_DIR/rsr/*.txt
ln -s -t $SDPC_INR_RUN_DIR/Staging/Left $INR_GOES_REFDATA_DIR/Reference/Left/*.nc
ln -s -t $SDPC_INR_RUN_DIR/Staging/Right $INR_GOES_REFDATA_DIR/Reference/Right/*.nc
#ln -s -t $SDPC_INR_RUN_DIR/CloudProducts $CLOUD_DATA_DIR/*.nc
