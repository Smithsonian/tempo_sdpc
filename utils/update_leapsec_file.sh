#! /bin/sh

set -u
set -e

LEAP_FILE="leap-seconds.list"

target=$(mktemp)
wget -q https://www.ietf.org/timezones/data/leap-seconds.list -O $target

if test -f $LEAP_FILE ; then
   old_md5sum=$(md5sum $LEAP_FILE | cut -d' ' -f1)
   new_md5sum=$(md5sum $target | cut -d' ' -f1)
   if test x"$old_md5sum" == x"$new_md5sum" ; then
      /bin/rm $target
      printf "update not needed\n"
      exit 0
  fi
  bkpfile="${LEAP_FILE}.$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  /bin/cp -a $LEAP_FILE $bkpfile
  printf "existing file copied to $bkpfile\n"
fi
chmod ugo+r $target
/bin/mv $target $LEAP_FILE
printf "updated $LEAP_FILE\n"
