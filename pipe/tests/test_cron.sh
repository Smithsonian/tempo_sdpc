#! /bin/sh

INCOMING="incoming"

if test -d "$INCOMING" ; then
  echo "directory exists: $INCOMING"
  exit 1
fi

. ./functions.sh

# Start the cache monitor.
# The incoming directory need not exist.
export SLANG_MODULE_PATH="../src"
../src/cachemon.sl ./cron.cfg -- > log.cron 2>&1 &
CACHEMON_PID="$!"

# Create the incoming directory and feed files
# matching the expected pattern.
/bin/mkdir "$INCOMING" || exit 1

tags="a b c"

for t1 in $tags ; do
 touch "$INCOMING/${t1}.tar"
done

wait_for_empty_dir $INCOMING 5

# Processing should be done by now.
kill -TERM $CACHEMON_PID
wait $CACHEMON_PID

# incoming directory should be empty at this point.
/bin/rmdir "$INCOMING" || exit 1

# cleanup
/bin/rm log.cron || exit 1
