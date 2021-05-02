#! /bin/sh

: "${SDPC_INR_RUN_DIR:?SDPC_INR_RUN_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

if test -e $SDPC_INR_RUN_DIR ; then
   echo "File exists: $SDPC_INR_RUN_DIR"
   exit 1
fi

INR_DIRS="config logs scantailoring \
          CloudProducts \
          Staging/Left \
          Staging/Right \
          Staging/rsr"

mkdir -p $SDPC_INR_RUN_DIR

for d in $INR_DIRS ; do
  mkdir -p $SDPC_INR_RUN_DIR/$d
done

ln -s $SDPC_RUN_DIR_MASTER/stage/granules/inr_input $SDPC_INR_RUN_DIR/Staging/Granules
ln -s $SDPC_RUN_DIR_MASTER/stage/granules/inr_output $SDPC_INR_RUN_DIR/Output
