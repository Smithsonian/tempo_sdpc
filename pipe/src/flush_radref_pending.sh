#! /bin/sh

if test -z "${SDPC_ROOT}"
then
  echo "The SDPC environment is not set"
  exit 1
fi

usage_message()
{
   echo "Usage: $(basename $0) [FILES]"
   exit 0
}

radref_search=$(config_setting radref.search)

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
   if test -f $path ; then
      # make sure radref_file is defined before we push this for processing
      if test -z "$radref_file" ; then
         printf "radref_file=\"$path\"\n" >> $tar_notice_path
      fi
      /bin/mv $tar_notice_path $SDPC_PIPE_DIR/stage/granules/level2_input
   fi
}

main()
{
  if test $# -ne 0 ; then
     case "$1" in
        --help) usage_message
          ;;
     esac
  fi

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

