#! /bin/sh

. ./functions.sh

if test -d incoming ; then
  echo "directory exists: incoming"
  exit 1
fi

if test -d logs ; then
  echo "directory exists: logs"
  exit 1
fi

if test -d out_multi ; then
  echo "directory exists: out_multi"
  exit 1
fi

export SLANG_MODULE_PATH=".."

/bin/mkdir logs || exit 1

NUM_CACHEMON=20
NUM_SCANS=50
NUM_GRANULES=10

# Start the cache monitors.
# The incoming directory need not exist.
pid_list=""
k=0
while test "$k" -lt $NUM_CACHEMON ; do
  ../cachemon.sl ./multi.cfg > logs/log.$k 2>&1 &
  pid_list="$! $pid_list"
  k=$(($k + 1))
done

# Create the incoming directory and feed files
# matching the expected pattern.
/bin/mkdir incoming || exit 1

s=0
while test $s -lt $NUM_SCANS ; do
   s3=$(printf "%03d" $s)
   g=0
    while test $g -lt $NUM_GRANULES ; do
      g3=$(printf "%03d" $g)
      touch "incoming/xxx_${s3}_${g3}.dat"
      g=$(($g + 1))
    done
  s=$(($s + 1))
done

wait_for_empty_dir incoming 5

# Processing should be done by now.
for pid in $pid_list ; do
  kill -0 $pid
  if test $? -eq 0 ; then
     kill -HUP $pid
     wait $pid
  else
     echo "cachemon pid=$pid already exited??"
  fi
done

# incoming directory should be empty at this point.
/bin/rmdir incoming || exit 1

# Count the files in each output directory
out_dirs=`ls out_multi`
for d in $out_dirs ; do
   num_files=`ls -1 out_multi/$d | wc -l`
   if ! test $num_files -eq $NUM_GRANULES ; then
      echo "*** ERROR: $d has $num_files files, expected $NUM_GRANULES"
      exit 1
   fi
done

# cleanup:
/bin/rm -rf out_multi || exit 1
/bin/rm -rf logs || exit 1
