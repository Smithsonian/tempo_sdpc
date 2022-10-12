#! /bin/sh

PGMNAME="$(basename $0)"

error_exit()
{
  echo "$PGMNAME [$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
  exit 1
}

usage_message()
{
   echo "Usage: $PGMNAME [options] <schedule-config>"
   echo "   Options:"
   echo "     --help        Print this usage message"
   echo "     --cfg         Print scan configuration template"
   echo "     --sch [args]  Print schedule configuration template."
   echo "                   Optional arguments: <maneuver-file> [<tailoring-file>]"
   echo "     --notes FILE  Text file containing special instructions."
   echo "                   (For details, see IOC/SDPC ICD TEMPO-09-0010)"
   exit 0
}

if test $# -eq 0 ; then
   usage_message
fi

# SDPC_OTS_ROOT is used in the plan.cfg file
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT is not set}"
: "${SDPC_ROOT:?SDPC_ROOT is not set}"

default_config_file="$SDPC_ROOT/share/plan.cfg"
example_sched_file="$SDPC_ROOT/share/plan_sched.sh.example"
notes_file=""

output_sched_file()
{
   maneuver_file="$1"
   shift
   if test $# -eq 1 ; then
      tailoring_file="$1"
   fi

   if ! test -f $example_sched_file ; then
      echo "*** Internal Error: no template: $example_sched_file"
      exit 1
   fi

   temp_file="$(mktemp)"
   /bin/cp $example_sched_file $temp_file

   if test -f "$maneuver_file" ; then
      tman_beg="$(grep table_begin_time $maneuver_file | cut -d, -f 1 | cut -d= -f 2)"
      tman_end="$(grep table_end_time $maneuver_file | cut -d, -f 1 | cut -d= -f 2)"
      sed -i -e s,"@MANEUVER_FILE_BEGIN@","$(date -u --date @$tman_beg +%Y-%m-%d)", \
             -e s,"@MANEUVER_FILE_END@","$(date -u --date @$tman_end +%Y-%m-%d)", \
             -e s,'maneuver_file=.*',"maneuver_file=\"$(realpath $maneuver_file)\"", \
             $temp_file
   elif ! test -z "$maneuver_file" ; then
      echo "*** Cannot access maneuver file: $maneuver_file"
      exit 1
   fi

   if test -f "$tailoring_file" ; then
      sed -i -e s,"tailoring_file=.*","tailoring_file=\"$tailoring_file\"", $temp_file
   elif ! test -z "$tailoring_file" ; then
      echo "*** Cannot access scan tailoring file: $tailoring_file"
      exit 1
   fi

   sed -i -e s,"@SCAN_START_DAY@","$(date --date friday +%Y-%m-%d)", \
          -e s,"@TAILOR_START_DAY@","$(date --date thursday +%Y-%m-%d)", \
          $temp_file

   cat $temp_file
   /bin/rm -f $temp_file
}

case "$1" in
   --sch)
        shift
        output_sched_file "$@"
        exit 0
        ;;

   --cfg)
        if test -f $default_config_file ; then
           cat $default_config_file
           exit 0
        fi
        ;;

   --help)
        usage_message
        ;;

   --notes)
        shift
        notes_file="$1"
        if ! test -r "$notes_file" ; then
           error_exit "*** Error: cannot read notes file: $notes_file"
        fi
        shift
        ;;

   *)
     # FALLTHRU
     ;;
esac

sched_file="$1"

if test -r "$sched_file" ; then
   . $(realpath $sched_file)
else
   error_exit "*** Error: cannot read schedule file: $sched_file"
fi

if test x"$config_file" = x ; then
   config_file="$default_config_file"
fi

: "${maneuver_file:?maneuver_file is not defined}"
: "${plan_type:?plan_type is not defined}"
: "${plan_num_days:?plan_num_days is not defined}"
: "${plan_start_day:?plan_start_day is not defined}"

trap error_exit ERR

# Get current time_t and convert to a UTC timestamp
timet=$(date +%s)
utc=$(date -u --date=@${timet} +%Y%m%dT%H%M%SZ)

plan_dirname="tempo_plan_${utc}"
target_dir="$plan_dirname"

if ! test -d $target_dir ; then
   mkdir -p $target_dir || error_exit "mkdir failed"
fi

echo "Writing output to $target_dir"

_tailor="${target_dir}/${utc}_scantailor.csv"
_master="${target_dir}/${utc}_masterscan.csv"
_plan="${target_dir}/${utc}_earthscan.csv"
_notes="${target_dir}/NOTES.txt"

# generate master scan table, and scan plan
plan -c $config_file -M $maneuver_file \
     -s $plan_type -d $plan_start_day -n $plan_num_days \
     $plan_options -o $_plan \
     -m $_master || error_exit "failed generating scan plan"

/bin/cp $config_file $target_dir/${utc}_plan.cfg
if ! test x"$notes_file" = x ; then
   /bin/cp $notes_file $_notes
else
   touch $_notes
fi

/bin/cp $maneuver_file $target_dir
_maneuver="$target_dir/$(basename $maneuver_file)"

if test x"$tailoring_file" = x ; then
   echo "WARNING: Creating synthetic scan tailoring file with zero offsets."
   if test \( x"${tailor_epoch}" = x \) -o \( x"${tailor_start_day}" = x \) -o \( x"${tailor_num_days}" = x \) ; then
      error_exit "Parameters for dummy scan tailoring file not specified"
   else
      plan -c $config_file --epoch $tailor_epoch \
        -d $tailor_start_day -n $tailor_num_days \
        -T $_tailor || error_exit "failed generating scan tailoring file"
   fi
elif test -r "$tailoring_file" ; then
  /bin/cp $tailoring_file $_tailor
else
  error_exit "cannot read $tailoring_file"
fi

# make generic local symlinks
ln -r -s $_plan $target_dir/earthscan.csv
ln -r -s $_master $target_dir/masterscan.csv
ln -r -s $_tailor $target_dir/scantailor.csv
ln -r -s $_maneuver $target_dir/maneuver.csv

tarfile="${plan_dirname}.tar.gz"

if test -d $target_dir ; then
  tar czf $tarfile -C $(dirname $target_dir) ${plan_dirname}
fi

echo "Created $tarfile"
