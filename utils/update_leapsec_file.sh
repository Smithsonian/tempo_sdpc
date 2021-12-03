#! /bin/sh

set -u
set -e

LEAP_FILE="leap-seconds.list"

target=$(mktemp)
wget -q https://www.ietf.org/timezones/data/leap-seconds.list -O $target

if test -f $LEAP_FILE ; then
  bkpfile="${LEAP_FILE}.$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  /bin/cp -a $LEAP_FILE $bkpfile
  printf "existing file copied to $bkpfile\n"
fi
chmod ugo+r $target
/bin/mv $target $LEAP_FILE
printf "updated $LEAP_FILE\n"
