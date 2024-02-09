#! /bin/sh

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

exit_usage()
{
   echo "Usage: $(basename $0) [options] <dirlist-file>"
   echo "   Options:"
   echo "   --help              Print this listing"
   exit "$1"
}

# Seems likely that the Level 0 version number won't need to change,
# but provide a mechanism anyway.
: "${SDPC_LEVEL0_VERSION:=1}"

run_l0_format()
{
   cachedir_list="$1"

   logdir="$SDPC_PIPE_DIR/log/level0"
   if ! test -d $logdir ; then
      echo "*** Error: cannot access log directory: $logdir"
      return
   fi
   echo "Running L0_format, log directory: $logdir"

   L0_format --archive "$SDPC_ARCHIVE_DIR" --register --verbose \
             --cache @${cachedir_list} \
             --Version "$SDPC_LEVEL0_VERSION" \
             --logdir "$logdir" \
             "$SDPC_PIPE_DIR/etc/l0_format.cfg" > "$logdir/l0_format.$$.log" 2>&1
}

main()
{
   if test $# -eq 0 || test $# -gt 1 || test x"$1" == x"--help" ; then
      exit_usage 0
   fi

   cachedir_list="$1"

   if ! test -f $cachedir_list ; then
      echo "*** Error: cannot read file: $cachedir_list"
      exit 1
   fi

   run_l0_format $cachedir_list
}

main "$@"
