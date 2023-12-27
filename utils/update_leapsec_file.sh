#! /bin/sh

set -u
set -e

LEAP_FILE="leap-seconds.list"

# 2023 Dec: Apparently only source for this file still works.
# Presumably this won't matter because in Nov 2022, BIPM resolved
# to phase out leap seconds by 2035. IERS prediction center says:
# "NO leap second will be introduced at the end of December 2023.
# The last leap second was positive and WAS introduced in UTC at the end of December 2016."

#URL="https://www.ietf.org/timezones/data/leap-seconds.list"
URL="https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list"

target=$(mktemp)
wget -q $URL -O $target

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
