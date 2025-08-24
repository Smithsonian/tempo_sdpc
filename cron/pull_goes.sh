#! /usr/bin/env bash

#set -e
set -u

if test $# -ne 2 ; then
    echo "Usage: $(basename $0) destdir USER@HOST"
    exit 0
fi

rootdir=$1
user_at_host=$2

# To simply tracking what we've downloaded,
# we first download files to an 'incoming' directory.
incoming_dir="$rootdir/incoming"
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir
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

# Move files to their final location, and
# register them in the asdc upload database

goes_archive.py --root $rootdir $incoming_dir/east_cmi
badchksum_east=$?

goes_archive.py --root $rootdir $incoming_dir/west_cmi
badchksum_west=$?

# Exit non-zero if any checksums failed
if test $badchksum_east -ne 0 -o $badchksum_west -ne 0 ; then
   exit 1
fi
