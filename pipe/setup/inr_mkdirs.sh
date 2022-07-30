#! /bin/sh

: "${SDPC_PIPE_DIR:?SDPC_PIPE_DIR not set, run this with sdpcrun.sh}"

set -u
set -e

inr_run_dir="$SDPC_PIPE_DIR/inr"

if test -e $inr_run_dir ; then
   echo "File exists: $inr_run_dir"
   exit 1
fi

INR_DIRS="config logs scantailoring CloudProducts Staging"

# GOES data will be delivered in these:
#    Staging/Left
#    Staging/Right
# which will be defined by symbolic links.
# During operations, these links will be updated daily

mkdir -p $inr_run_dir

for d in $INR_DIRS ; do
  mkdir -p $inr_run_dir/$d
done

ln -s $SDPC_PIPE_DIR/stage/granules/inr_input $inr_run_dir/Staging/Granules
ln -s $SDPC_PIPE_DIR/stage/granules/inr_output $inr_run_dir/Output
