#! /bin/bash

set -u
set -e

if test $# -eq 0 ; then
   echo "Usage:  $0 <filedb-executable-path>"
   exit 1
fi

FILEDB_EXEC="$1"

TEST_DBNAME="met:nam227"
SEARCH_TIME="1478107809"   # UTC sec since unix epoch

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
  __tic=$(date +"%s%N")
}

toc()
{
  __toc=$(date +"%s%N")
}

assert_duration_exceeds()
{
   expect=$1

   duration=$(echo "scale=9; ($__toc - $__tic)/10^9" | bc)
   value=$(echo $duration \> $expect | bc)

   if test $value -eq 0 ; then
      printf "*** Duration %g sec is less than expected value %g sec\n" $duration $expect
      exit 1
   fi
}

assert_duration_lessthan()
{
   expect=$1

   duration=$(echo "scale=9; ($__toc - $__tic)/10^9" | bc)
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

# Since leap seconds aren't critical, we search on UTC time stamps based on the Unix epoch.
# Defining search_time = time since the unix epoch 1970-01-01T00:00:00 UTC we have:
# 1478111409  /tmp/filedb_test.houck/ancillary/met/nam227/2016/11/02/2016110218.nam.tffz.conusnest.hiresf18.tm00.grib2
# 1478179809  /tmp/filedb_test.houck/ancillary/snow/nsidc/2016/11/NISE_SSMISF18_20161103.HDFEOS
# 1478136609  /tmp/filedb_test.houck/archive/L1/irr/2016/11/03/TEMPO_irr_L1_V01_20161103T013009Z.nc
# 1478086209  /tmp/filedb_test.houck/archive/L0/drk/2016/11/02/TEMPO_drk_L0_V01_20161102T113009Z.nc

check_search_result met:nam227 1478111409
check_search_result snow:nsidc 1478179809
check_search_result tempo:irr  1478136609
check_search_result tempo:drk  1478086209

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
