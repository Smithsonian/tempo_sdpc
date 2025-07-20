#! /bin/sh

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

if test $# -eq 0 ; then
   echo "Usage: $(basename $0) [radref_search]"
   echo " Where radref_search = 0|1"
   exit 0
fi

radref_search=$1
case "$radref_search" in
  0 | 1 )
        # valid input
        ;;
  * )
    echo "*** Error: invalid input: $radref_search (input value must be 0 or 1)"
    exit 1
    ;;
esac

process_tar_file()
{
   tar_notice_path="$1"

   radref_file=""
   . $(realpath $tar_notice_path)

   if test $radref_search -ne 0 ; then
      path=$(select_radref.py $rad_filename)
   else
      path=$(select_radref.py --thisscan $rad_filename)
   fi
   if test -n "$path" && test -f "$path" ; then
      # make sure radref_file is defined before we push this for processing
      if test -z "$radref_file" ; then
         printf "radref_file=\"$path\"\n" >> $tar_notice_path
      fi
      /bin/mv $tar_notice_path $SDPC_PIPE_DIR/stage/granules/level2_input
   fi
}

main()
{
  pending_dir="$SDPC_PIPE_DIR/stage/granules/level2_input/radref_pending"
  if ! test -d $pending_dir ; then
     exit 0
  fi

  files=$(find $pending_dir -mindepth 1 -maxdepth 1 -type f)

  for f in $files ; do
     process_tar_file $f
  done
}

main "$@"

