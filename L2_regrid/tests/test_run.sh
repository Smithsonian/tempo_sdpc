#! /bin/sh

# Test data and L2_regrid output goes into this directory
TMPDIR="tmp"

# Note that test_run.cfg should be consistent with these values
NUM_STEPS=40
NUM_PIXELS=100

mkdir -p $TMPDIR
./mktestdata_l2 -o $TMPDIR -m $NUM_STEPS -p $NUM_PIXELS
../src/gnuobjs/L2_regrid test_run.cfg
