#! /bin/sh

set -u
set -e

LEAP_FILE="leap-seconds.list"

if test -f $LEAP_FILE ; then
  suffix=$(date --iso-8601)
  newfile="${LEAP_FILE}.${suffix}"
  /bin/mv $LEAP_FILE $newfile
  printf "moved existing file to $newfile\n"
fi

wget -q https://www.ietf.org/timezones/data/leap-seconds.list
