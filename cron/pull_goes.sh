#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

if ! command -v sha1sum &> /dev/null
then
    echo "ERROR: sha1sum could not be found"
    exit 1
fi

#set -e
set -u

if test $# -ne 1 ; then
    echo "Usage: $0 USER@HOST"
    exit 0
fi
user_at_host=$1

rootdir="${SDPC_ANCILLARY_ROOT}/var/goes"

# It's convenient if a single directory contains
# all GOES imagery needed for a single operational day.
# For this reason, we compute the day-of-year using
# the satellite-local time zone for TEMPO.
# Weirdly, POSIX requires a positive TZ offset
# for time zones *west* of the prime meridian.
subdir="$(TZ='UTC+6' date +%Y/%j)"

# To simply tracking what we've downloaded,
# we first download files to an 'incoming' directory.
incoming_dir="$rootdir/incoming"
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir
fi

# After database entry the files will be moved
# to their final location:
target_dir="$rootdir/$subdir"
if ! test -d $target_dir ; then
   mkdir -p $target_dir
fi

lftp $user_at_host <<- EOF
   set log:file/xfer ""
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set ssl:verify-certificate no
   set mirror:require-source true
   set mirror:sort-by name
   set mirror:order *.sha1 *.nc
   mirror -c -O $incoming_dir -F /PDAFileLinks/??st_cmi
   quit
EOF

verify_cksum()
{
   nc_path=$1
   expected_cksum=$(cat ${nc_path}.sha1)
   actual_cksum=$(sha1sum $nc_path | cut -d' ' -f1)
   if test x"$expected_cksum" = x"$actual_cksum" ; then
      #echo "$(date -u +%Y%m%d%H%M%SZ): good checksum: $nc_path"
      return 0
   else
      echo "$(date -u +%Y%m%d%H%M%SZ): WARNING: bad checksum: $nc_path"
      return 1
   fi
}
export -f verify_cksum

# Move files to their final location, and
# register them in the asdc upload database

move_and_register_files()
{
   dir="$1"
   dbfile="$2"
   files=$(find $incoming_dir/$dir -name "*.nc")

   dest_dir="$target_dir/$dir"
   # make sure the destination dir exists
   if ! test -d $dest_dir ; then
       mkdir -p $dest_dir
   fi

   # files with bad checksums go here
   badchksum_dir="$dest_dir/bad_checksum"
   badchksum=0

   for f in $files ; do
      if test -f ${f}.sha1 ; then
         # verify checksum
         verify_cksum $f
         bad=$?
         if test $bad -ne 0 ; then
            badchksum=1
            mkdir -p $badchksum_dir
            /bin/mv $f ${f}.sha1 $badchksum_dir
            continue
         fi
      fi
      bn=$(basename $f)
      final_path="$dest_dir/$bn"
      /bin/mv $f $final_path
      asdc_files.py --dbfile $dbfile --add $final_path
      # move the sha1 file too
      if test -f ${f}.sha1 ; then
         /bin/mv ${f}.sha1 $dest_dir
      fi
   done

   return $badchksum
}

move_and_register_files east_cmi $rootdir/cmieast.sqlite
badchksum_east=$?

move_and_register_files west_cmi $rootdir/cmiwest.sqlite
badchksum_west=$?

# The 'today' symlink always points to today's GOES data
cd $rootdir
ln -nf -s $subdir today

# Exit non-zero if any checksums failed
if test $badchksum_east -ne 0 -o $badchksum_west -ne 0 ; then
   exit 1
fi
