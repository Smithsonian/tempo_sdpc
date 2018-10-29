#! /bin/bash

set -u
set -e

if test $# -eq 0 ; then
   echo "Usage:  $0 <filedb-executable-path>"
   exit 1
fi

FILEDB_EXEC="$1"

TEST_DBNAME="met:nam227"
SEARCH_TIME="541647009"

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
      printf "Duration %g sec is less than expected value %g sec\n" $duration $expect
      exit 1
   fi
}

assert_duration_lessthan()
{
   expect=$1

   duration=$(echo $__toc - $__tic | bc)
   value=$(echo $duration \< $expect | bc)

   if test $value -eq 0 ; then
      printf "Duration %g sec is longer than expected value %g sec\n" $duration $expect
      exit 1
   fi
}

#-------------- Exercise database options

DBNAME_LIST="met:nam227 snow:nsidc tempo:irr tempo:drk"
for db in $DBNAME_LIST ; do
   $FILEDB_EXEC $db --update
done
for db in $DBNAME_LIST ; do
   $FILEDB_EXEC $db --find --sec $SEARCH_TIME > /dev/null
done

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
