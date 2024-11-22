#! /bin/sh

check_cmd_exists()
{
   cmd="$1"
   if ! command -v $cmd &> /dev/null ; then
      echo "*** Error: $cmd not found"
      exit 1
   fi
}

check_cmd_exists multitail
check_cmd_exists whiptail

exit_usage()
{
   echo "Usage: $(basename $0) [options] [service [service] ...]"
   echo "   Options:"
   echo "   --help     Print this listing"
   echo "   --dir      Primary node pipeline directory [default=\$SDPC_PIPE_DIR]"
   echo "   --ioc      Show IOC interface services"
   echo "   --sci      Show science processing services"
   echo "   --nrt      Show NRT services"
   echo "   --menu     Display service menu"
   exit "$1"
}

default_logdirs()
{
  logdirs="level1a inr level1b level2 level3 register asdc"
}

choose_logdirs()
{
   read -r -d '' itemlist <<-'EOF'
	iocpull off
	iocpullraw off
	level0 off
	level1a on
	inr on
	level1b on
        level1b_nrt1 off
        level1b_nrt2 off
	level2 on
        level2_nrt off
	level3 on
	trend off
	register on
	asdc on
	pipecron off
	EOF

  num_items=$(echo $itemlist | wc -w)
  num_items=$((num_items / 2))

  num_lines=$((num_items + 10))
  num_cols=32

  logdirs=$(
  whiptail --checklist "Choose which logs to view:" \
           --title "SDPC services" \
           --noitem \
           $num_lines $num_cols $num_items $itemlist 3>&2 2>&1 1>&3
  )
  status="$?"
  if test $status = 1; then
     exit 0
  elif ! test $status = 0 ; then
     exit $status
  fi
}

main()
{
  default_logdirs

  if ! test -z "$SDPC_PIPE_DIR" ; then
     root="$SDPC_PIPE_DIR"
  fi

  # Process optional args
  while [ "$#" != "0" ]
  do
     case "$1" in
       --*)
         case "$1" in
           --help)
             exit_usage 0
             ;;
	   --dir)
	      shift
 	      root="$1"
	      shift
              if ! test -d "$root" ; then
                 echo "Nonexistent directory: $root"
                 exit 1
              fi
	      ;;
           --ioc)
             shift
             logdirs="iocpull iocpullraw level0 register asdc"
             ;;
           --nrt)
             shift
             logdirs="level1b_nrt1 level1b_nrt2 level2_nrt register asdc"
             ;;
           --sci | --science)
             shift
             logdirs="level1a inr level1b level2 level3 trend register asdc"
             ;;
           --menu)
              shift;
              choose_logdirs
              ;;
           *)
             echo "Unknown option: $1"
             exit 1
             ;;
         esac
         ;;
       *)
         logdirs="$@"
         break
         ;;
     esac
  done

  logfiles=$(printf "$root/log/%s/current " $logdirs | tr -d \")

  multitail -o default_convert:qmailtimestr:current \
          -v -E TIO_TRACE $logfiles
}

main "$@"
