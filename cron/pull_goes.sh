#! /usr/bin/env bash

: "${SDPC_ANCILLARY_ROOT:?SDPC_ANCILLARY_ROOT not set}"

set -e
set -u

rootdir="${SDPC_ANCILLARY_ROOT}/goes"

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

lftp AO_TEMPO_SERV1@140.90.190.143 <<- EOF
   set xfer:log-file
   set xfer:use-temp-file yes
   set xfer:temp-file-name *.lftp
   set ssl:verify-certificate no
   set mirror:require-source true
   set mirror:sort-by name
   set mirror:order *.sha1 *.nc
   mirror -c -O $incoming_dir -F /PDAFileLinks/g1?_cmi
   quit
EOF

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

   for f in $files ; do
      bn=$(basename $f)
      final_path="$dest_dir/$bn"
      /bin/mv $f $final_path
      bin/asdc_files.py --dbfile $dbfile --add $final_path
      # move the sha1 files too
      if test -f ${f}.sha1 ; then
         /bin/mv ${f}.sha1 $dest_dir
      fi
   done
}

move_and_register_files g16_cmi $rootdir/cmig16.sqlite
move_and_register_files g17_cmi $rootdir/cmig17.sqlite

# The 'today' symlink always points to today's GOES data
cd $rootdir
ln -nf -s $subdir today
