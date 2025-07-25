# This file is intended to contain only function definitions,
# so running it alone doesn't do anything.
# The functions are intended to be called from a parent script
# which imports these definitions.

PROGNAME="$(basename $0)"
mdc_error_exit()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error: line $2: $1"
    exit 1
  fi
}
trap 'mdc_error_exit $? $LINENO' EXIT

make_path_for_cldo4_destripe_file()
{
   first_path="$1"

   destripe_filename="$(basename $first_path | sed -E -e s,CLDO4,DSTRCLDO4, -e 's,T[0-9]+Z_S[0-9]+G[0-9]+,,')"

   day_dir=$(level1_info -l $first_path)
   destripe_dir="$SDPC_ARCHIVE_DIR/L2/RAD/$day_dir/destripe/CLDO4"

   echo "$destripe_dir/$destripe_filename"
}

make_cldo4_destripe_file()
{
   l2_path_list_file="$1"

   first_path=$(head -1 $l2_path_list_file)

   destripe_path="$(make_path_for_cldo4_destripe_file $first_path)"
   destripe_dir=$(dirname $destripe_path)
   if ! test -d $destripe_dir ; then
      mkdir -p $destripe_dir
   fi

   log_file="$destripe_dir/make_destripe.CLDO4.log"

   # Generate the destriping correction file
   tempo_destripe_regular.py --mode derive --desfnm $destripe_path --list4descor $l2_path_list_file > $log_file 2>&1
   if test $? -ne 0 ; then
      mdc_error_exit "make_destripe.py failed (see $log_file)" $LINENO
   fi

   # Register the file in the sqlite database.
   ln -s $destripe_path $SDPC_ARCHIVE_DIR/registry/incoming
}
