#! /bin/sh

# Test data and L2_regrid output goes into this directory
TMPDIR="tmp"
L2_REGRID="$1"
PREFIX=""

# Note that test_run.cfg should be consistent with these values
NUM_STEPS=40
NUM_PIXELS=100

mkdir -p $TMPDIR
./mktestdata_l2 -o $TMPDIR -m $NUM_STEPS -p $NUM_PIXELS
$PREFIX $L2_REGRID --ignore test_run.cfg
