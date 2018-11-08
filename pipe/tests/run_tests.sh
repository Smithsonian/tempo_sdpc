#! /bin/sh

TESTS="$@"

for t in $TESTS ; do
  echo "Running $t:"
  ./$t
  if test "$?" -ne 0 ; then
     echo "FAILED"
     exit 1
  fi
  echo "ok"
done

echo "All tests passed"
