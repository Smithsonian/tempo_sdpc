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

process_tar_file()
{
   tar_notice_path="$1"

   . $(realpath $tar_notice_path)
   path=$(select_radref.py --thisscan $rad_filename)
   if test -f $path ; then
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

  files=$(find $SDPC_PIPE_DIR/stage/granules/level2_input/radref_pending -mindepth 1 -maxdepth 1 -type f)

  for f in $files ; do
     process_tar_file $f
  done
}

main "$@"

