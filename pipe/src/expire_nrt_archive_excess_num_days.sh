#! /bin/sh

# Delete only directories with names that match this pattern
_pattern="D?????"

# Delete things only in directories that contain a file with this name:
_enable_filename="AUTO_DELETE_ENABLED"

_pgm_name="$(basename $0)"

usage()
{
   echo "Usage: $_pgm_name num_keep"
   echo "       Note: num_keep > 1 is required."
   exit $1
}

if test $# -ne 1 ; then
   usage 0
fi

: "${SDPC_ARCHIVE_DIR:?SDPC_ARCHIVE_DIR not set}"

_num_keep=$1
shift

# ensure _num_keep is an integer > 1
if test "$_num_keep" -le 1 ; then
   usage 1
fi

expire_pattern_dirs()
{
   rootdir=$1

   # Reject relative paths
   if test x"$rootdir" != x"${rootdir#*..}" ; then
      echo "WARNING: directory path contains '..': $rootdir"
      return
   fi

   if ! test -d $rootdir ; then
      return
   fi

   enable_file="$rootdir/$_enable_filename"

   # Operate only in specially marked directories
   if test -f "$enable_file" ; then
      # By default, find does not follow symlinks. That's the safest choice here.
      find $rootdir -mindepth 1 -maxdepth 1 -type d -name "$_pattern" | sort | head -n -$_num_keep | xargs --no-run-if-empty /bin/rm -rf
   else
      echo "WARNING: directory expiration not enabled: $rootdir (to enable, touch $enable_file)"
   fi
}

main()
{
   # NRT archive directory tree must exist
   nrt_archdir="$SDPC_ARCHIVE_DIR/NRT"
   if ! test -d "$nrt_archdir" ; then
      echo "WARNING:  nonexistent directory: $nrt_archdir"
      exit 0
   fi

   # Only operate on NRT archive subdirectories
   subdirs="L1/RAD L2/RAD L3/RAD"

   for d in $subdirs ; do
       expire_pattern_dirs "$nrt_archdir/$d"
   done
}

main "$@"
