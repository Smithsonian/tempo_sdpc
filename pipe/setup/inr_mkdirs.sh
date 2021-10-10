#! /bin/sh

: "${SDPC_RUN_DIR_INR:?SDPC_RUN_DIR_INR not set, run this with sdpcrun.sh}"

set -u
set -e

if test -e $SDPC_RUN_DIR_INR ; then
   echo "File exists: $SDPC_RUN_DIR_INR"
   exit 1
fi

INR_DIRS="config logs scantailoring CloudProducts Staging/rsr"

# GOES data will be delivered in these:
#    Staging/Left
#    Staging/Right
# which will be defined by symbolic links.
# During operations, these links will be updated daily

mkdir -p $SDPC_RUN_DIR_INR

for d in $INR_DIRS ; do
  mkdir -p $SDPC_RUN_DIR_INR/$d
done

ln -s $SDPC_RUN_DIR_MASTER/stage/granules/inr_input $SDPC_RUN_DIR_INR/Staging/Granules
ln -s $SDPC_RUN_DIR_MASTER/stage/granules/inr_output $SDPC_RUN_DIR_INR/Output
