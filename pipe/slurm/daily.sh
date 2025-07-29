#! /bin/sh

# 0. This script is run by cachemon, upon completion of all level 2
#    data products of a particular type from a single day.
#    The script is triggered by the appearance of a file in
#      $SDPC_PIPE_DIR/stage/days
#    The argument passed to the script is the path to this file.

set -u

if test -z SDPC_ROOT ; then
   printf "SDPC environment is not set\n"
   exit 1
fi

if test $# -ne 1 ; then
   echo "Usage: $0 <L2-pathlist-file>"
   exit 1
fi

pathlist_file="$1"

PROGNAME="$(basename $0)"
catch()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error $1 occurred on $2"
  fi
}
trap 'catch $? $LINENO' EXIT

log_message()
{
   printf "${PROGNAME}: $1\n"
}

error_exit()
{
   printf "*** Error: ${PROGNAME}: $1\n"
   exit 1
}

pathlist_basename=$(basename $pathlist_file | sed -e s"/^[.]//")
log_message "processing $pathlist_basename"

# Import functions to generate destriping correction files
. $SDPC_ROOT/bin/make_destripe.sh
. $SDPC_ROOT/bin/make_destripe_cldo4.sh

# Parse pathlist filename:
pathlist_basename_sans_extname=$(basename $pathlist_file .lis)
product_type=$(echo $pathlist_basename_sans_extname | cut -d_ -f2)

# Make daily destriping files for selected products
if test -n $SDPC_MAKE_DESTRIPE_TG ; then
   _destripe_products=$(echo $SDPC_MAKE_DESTRIPE_TG | tr , ' ')
   for p in $_destripe_products ; do
       if test $p = $product_type ; then
          make_day_destripe_file $pathlist_file
          break
       fi
   done
fi

case "$product_type" in
   CLDO4_L2 )
   make_cldo4_destripe_file $pathlist_file
   ;;

   * )
   ;;
esac

/bin/rm -f "$pathlist_file"
