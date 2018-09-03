#! /bin/sh

INCOMING="incoming"

if test -d "$INCOMING" ; then
  echo "directory exists: $INCOMING"
  exit 1
fi

. ./functions.sh

# Start the cache monitors.
# The incoming directory need not exist.
export SLANG_MODULE_PATH=".."
../cachemon.sl ./stage1.cfg -- 3 > log.stage1 2>&1 &
CACHEMON_PID_STAGE1="$!"

../cachemon.sl ./stage2.cfg > log.stage2 2>&1 &
CACHEMON_PID_STAGE2="$!"

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
kill -HUP $CACHEMON_PID_STAGE1
wait $CACHEMON_PID_STAGE1

kill -HUP $CACHEMON_PID_STAGE2
wait $CACHEMON_PID_STAGE2

# incoming directory should be empty at this point.
/bin/rmdir "$INCOMING" || exit 1

/bin/ls -R out_stage2 > out_stage2.lst
lst2_got=`cat out_stage2.lst | sort`
lst2_ok=`cat expected-out_stage2.lst | sort`
if ! test X"$lst2_got" = X"$lst2_ok" ; then
   exit 1
fi

# cleanup:

# out_stage1 should be empty
/bin/rmdir out_stage1 || exit 1
/bin/rm -rf out_stage2 || exit 1
/bin/rm out_stage2.lst log.stage2 log.stage1 || exit 1
