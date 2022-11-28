#! /bin/sh

when="7 days ago"
when_timet=""
verbose=0

asdc_status_table()
{
   sqlite_file="$1"

   sql_newest="select timet,(now-timet)/3600.0 from (select max(strftime('%s',timestamp)) as timet,strftime('%s',current_timestamp) as now from File_Table)"

   if test $verbose -ne 0 ; then
      filter=""
   else
      filter="asdc_status != 3 and"
   fi

   # Note strftime()>t didn't work, but strftime()-t>0 did.  sqlite3 bug I guess.
   sql="select date(timestamp),asdc_status,count(all) from File_Table \
where $filter strftime('%s',timestamp) - $when_timet > 0 \
group by date(timestamp),asdc_status"

   if test -f $sqlite_file ; then
      newest_timet_duration=$(sqlite3 -separator , $sqlite_file "$sql_newest")
      result=$(sqlite3 -header $sqlite_file "$sql")
   else
      result="*** NOT FOUND ***"
      newest_timet_duration=""
   fi

   printf "\n dbfile: $sqlite_file\n"
   if test -n "$newest_timet_duration" ; then
      newest_timet=$(echo $newest_timet_duration | cut -d, -f1)
      duration_hours=$(echo $newest_timet_duration | cut -d, -f2)
      printf " newest: %0.2f hours old (%s)\n" "$duration_hours" "$(date -u --date @$newest_timet +%Y-%m-%dT%H:%M:%SZ)"
      printf "uploads:"
      if test -z "$result" ; then
        printf " OK"
      fi
      printf "\n"
   fi
   if test -n "$result" ; then
      printf "$result\n" | sed -e "s,^,\t,"
   fi
}

exit_usage()
{
  status=$1
  printf "Usage: %s [options]\n" $(basename $0)
  printf "Options:\n"
  printf "   --help        Print this usage message\n"
  printf "   --start WHEN  Start time to be interpreted by 'date'\n"
  printf "   -v            Generate verbose output\n"
  exit $status
}

main()
{
   # Process optional args
   while [ "$#" != "0" ]
   do
     case "$1" in
       -v)
        verbose=1
        shift
        ;;
       --*)
         case "$1" in
           --help)
             exit_usage 0
             ;;
           --start)
              shift;
              when="$1"
              shift;
              ;;
           *)
             echo "Unknown option: $1"
             exit 1
             ;;
         esac
         ;;
       *)
         break
         ;;
     esac
   done

   if test "$#" -ne 0 ; then
      echo "*** Error: unexpected command line arguments"
      echo "$@"
      exit_usage 1
   fi

   when_utc=$(date -u --date "$when" +%Y-%m-%dT%H:%M:%SZ)
   when_timet=$(date --date "$when" +%s)

   echo
   printf "ANCILLARY DATA STATUS\n"
   printf " start: $when_utc\n"
   printf "\n======= GOES Imagery:\n"
   asdc_status_table "$SDPC_ANCILLARY_ROOT/var/goes/cmieast.sqlite"
   asdc_status_table "$SDPC_ANCILLARY_ROOT/var/goes/cmiwest.sqlite"
   printf "\n======= GEOS Composition Forecasts:\n"
   asdc_status_table "$SDPC_ANCILLARY_ROOT/var/geoscf/geoscf.sqlite"
   printf "\n======= IMS Snow & Ice Cover Maps:\n"
   asdc_status_table "$SDPC_ANCILLARY_ROOT/var/ims/ims.sqlite"
}

main "$@"
