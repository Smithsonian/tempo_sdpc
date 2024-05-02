#! /bin/sh

PGMNAME="$(basename $0)"

# SDPC_OTS_ROOT is used in the plan.cfg file
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT is not set}"
: "${SDPC_ROOT:?SDPC_ROOT is not set}"

config_file="$SDPC_ROOT/share/plan.cfg"
points_file="$SDPC_ROOT/share/report_points.csv"
plan_pgm="$SDPC_ROOT/bin/plan"
epoch="1980-01-06T00:00:00Z"

error_exit()
{
  echo "$PGMNAME [$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
  exit 1
}
trap error_exit ERR

usage_message()
{
   echo "Usage: $PGMNAME [options] YYYY-MM-DD YYYY-MM-DD"
   echo "   Options:"
   echo "     --help          Print this usage message"
   echo "     --config FILE   Plan config file [default: $config_file]"
   echo "     --points FILE   Points file [default: $points_file]"
   exit 0
}

main ()
{
   if [ "$#" = "0" ]; then
     usage_message
   fi

   # Process optional args
   while [ "$#" != "0" ]
   do
     case "$1" in
       --*)
         case "$1" in
           --help)
             usage_message
             ;;
           --config)
              shift;
              config_file="$1"
              shift;
              ;;
           --points)
              shift;
              points_file="$1"
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

  if test "$#" -ne 2 ; then
    usage_message 1
  fi

  if ! test -r $config_file ; then
     echo "*** Error: cannot read config file: $config_file"
     exit 1
  fi
  if ! test -r $points_file ; then
     echo "*** Error: cannot read points file: $points_file"
     exit 1
  fi

  yyyy_mm_dd1="$1"
  yyyy_mm_dd2="$2"
  s0=$(date --date $yyyy_mm_dd1 +%s)
  s1=$(date --date $yyyy_mm_dd2 +%s)
  # shell only does integer arithmetic
  num_days=$((1 + ($s1 - $s0 + 43200)/86400))

  report_file="tempo_report_${yyyy_mm_dd1}_to_${yyyy_mm_dd2}.csv"
  tailor_file="tempo_tailor_${yyyy_mm_dd1}_to_${yyyy_mm_dd2}.csv"

  common_args="--config $config_file --epoch $epoch --date ${yyyy_mm_dd1} --ndays $num_days"

  echo "Generating $report_file"
  ($plan_pgm $common_args --points $points_file > $report_file) || error_exit "failed generating reporting file"

  echo "Generating $tailor_file"
  $plan_pgm $common_args -T $tailor_file || error_exit "failed generating scan tailoring file"
}

main "$@"

