#! /bin/sh

INCOMING="incoming"

if test -d "$INCOMING" ; then
  echo "directory exists: $INCOMING"
  exit 1
fi

. ./functions.sh

# Start the cache monitor.
# The incoming directory need not exist.
export SLANG_MODULE_PATH=".."
../cachemon.sl ./process1.cfg > log.process1 2>&1 &
CACHEMON_PID="$!"

# Create the incoming directory and feed files
# matching the expected pattern.
/bin/mkdir "$INCOMING" || exit 1

tags="111 222 333"

for t1 in $tags ; do
   for t2 in $tags ; do
   touch "$INCOMING/xxx_${t1}_${t2}.dat"
   done
done

wait_for_empty_dir $INCOMING 5

# Processing should be done by now.
kill -HUP $CACHEMON_PID
wait $CACHEMON_PID

# incoming directory should be empty at this point.
/bin/rmdir "$INCOMING" || exit 1

# Check the output directory contents:
/bin/ls -R out_process1 > out_process1.lst
lst_got=`cat out_process1.lst | sort`
lst_ok=`cat expected-out_process1.lst | sort`
if ! test X"$lst_got" = X"$lst_ok" ; then
   exit 1
fi

# cleanup
/bin/rm -rf out_process1 || exit 1
/bin/rm out_process1.lst log.process1 || exit 1
