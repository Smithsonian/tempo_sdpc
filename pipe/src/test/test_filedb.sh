#! /bin/bash

set -u
set -e

if test $# -eq 0 ; then
   echo "Usage:  $0 <filedb-executable-path>"
   exit 1
fi

FILEDB_EXEC="$1"

TEST_DBNAME="met:nam227"
SEARCH_TIME="531379809"

DELAY=5
APPROX_EXPECTED_DELAY_WHEN_BLOCKED=3

exit_message="FAIL";

print_exit_message()
{
   printf "$exit_message\n"
}
trap print_exit_message EXIT ERR

tic()
{
  __tic=$(date +"%s")
}

toc()
{
  __toc=$(date +"%s")
}

assert_duration_exceeds()
{
   expect=$1

   duration=$(echo $__toc - $__tic | bc)
   value=$(echo $duration \> $expect | bc)

   if test $value -eq 0 ; then
      printf "*** Duration %g sec is less than expected value %g sec\n" $duration $expect
      exit 1
   fi
}

assert_duration_lessthan()
{
   expect=$1

   duration=$(echo $__toc - $__tic | bc)
   value=$(echo $duration \< $expect | bc)

   if test $value -eq 0 ; then
      printf "*** Duration %g sec is longer than expected value %g sec\n" $duration $expect
      exit 1
   fi
}

check_search_result()
{
   db=$1
   expect_tt=$2

   fn=$($FILEDB_EXEC $db --find --sec $SEARCH_TIME)
   tt=$(cat $fn)
   #printf "$tt $fn\n"

   if test $tt -ne $expect_tt ; then
      printf "*** unexpected search result: %ld != %ld\n" $tt $expect_tt
      exit 1
   fi
}

printf "Test: database initialization:\n"

DBNAME_LIST="met:nam227 snow:nsidc tempo:irr tempo:drk"
for db in $DBNAME_LIST ; do
   $FILEDB_EXEC $db --update
done

printf "Test: database search result accuracy:\n"

# SEARCH_TIME="531379809"
# 531379809 /tmp/filedb_test.jhouck/ancillary/met/nam227/2016/11/02/2016110217.nam.tffz.conusnest.hiresf17.tm00.grib2
# 531365409 /tmp/filedb_test.jhouck/ancillary/snow/nsidc/2016/11/NISE_SSMISF18_20161102.HDFEOS
# 531322209 /tmp/filedb_test.jhouck/archive/L1/irr/2016/11/02/TEMPO_irr_L1_V01_20161102T013009Z.nc
# 531358209 /tmp/filedb_test.jhouck/archive/L0/drk/2016/11/02/TEMPO_drk_L0_V01_20161102T113009Z.nc

check_search_result met:nam227 531379809
check_search_result snow:nsidc 531365409
check_search_result tempo:irr  531322209
check_search_result tempo:drk  531358209

#-------------- Test database file locking

printf "Test: read should block until ongoing writes are completed:\n"
#-----------------------------------------------------------------------
# Initialize LUT, but delay lock release.
# While giving that process enough time to get the lock,
# try performing a lookup while the LUT is locked.
# The 2nd process should block until the first process completes
# and releases the lock:
$FILEDB_EXEC -d $DELAY $TEST_DBNAME --update &
sleep 0.1
tic
$FILEDB_EXEC $TEST_DBNAME --find --sec $SEARCH_TIME > /dev/null
toc
assert_duration_exceeds $APPROX_EXPECTED_DELAY_WHEN_BLOCKED

printf "Test: write should block until ongoing writes are completed:\n"
#-----------------------------------------------------------------------
# Initialize LUT, but delay lock release.
# While giving that process enough time to get the lock,
# try running a second initialization process, while the LUT is locked.
# The 2nd process should block until the first process completes
# and releases the lock:
$FILEDB_EXEC -d $DELAY $TEST_DBNAME --update &
sleep 0.1
tic
$FILEDB_EXEC $TEST_DBNAME --update
toc
assert_duration_exceeds $APPROX_EXPECTED_DELAY_WHEN_BLOCKED

printf "Test: write should block until ongoing reads are completed:\n"
#-----------------------------------------------------------------------
# Perform a lookup, but delay read-lock release.
# While giving that process enough time to get the lock,
# try initializing the LUT while the LUT is locked.
# The 2nd process should block until the first process completes
# and releases the lock:
$FILEDB_EXEC $TEST_DBNAME -d $DELAY --find --sec $SEARCH_TIME > /dev/null &
sleep 0.1
tic
$FILEDB_EXEC $TEST_DBNAME --update
toc
assert_duration_exceeds $APPROX_EXPECTED_DELAY_WHEN_BLOCKED

printf "Test: simultaneous reads are OK and do not block each other:\n"
#-----------------------------------------------------------------------
# Perform a lookup, but delay read-lock release.
# While giving that process enough time to get the lock,
# try performing another lookup while the LUT is locked.
# The 2nd process should not block.
$FILEDB_EXEC $TEST_DBNAME -d $DELAY --find --sec $SEARCH_TIME > /dev/null &
sleep 0.1
tic
$FILEDB_EXEC $TEST_DBNAME --find --sec $SEARCH_TIME > /dev/null
toc
assert_duration_lessthan 1

exit_message="OK"
