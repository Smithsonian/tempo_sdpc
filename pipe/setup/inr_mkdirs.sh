#! /bin/sh

: "${SDPC_INR_RUN_DIR:?SDPC_INR_RUN_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

if test -e $SDPC_INR_RUN_DIR ; then
   echo "File exists: $SDPC_INR_RUN_DIR"
   exit 1
fi

INR_DIRS="config logs scantailoring \
          Output CloudProducts \
          Staging/Granules \
          Staging/Left \
          Staging/Right \
          Staging/rsr"

mkdir -p $SDPC_INR_RUN_DIR

for d in $INR_DIRS ; do
  mkdir -p $SDPC_INR_RUN_DIR/$d
done
