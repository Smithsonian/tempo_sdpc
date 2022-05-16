#! /bin/sh

PGMNAME="$(basename $0)"

error_exit()
{
  echo "$PGMNAME [$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
  exit 1
}

usage_message()
{
   echo "Usage: $PGMNAME <schedule-config> | [options]"
   echo "   Options:"
   echo "     --help        Print this usage message"
   echo "     --sch         Print schedule configuration template"
   echo "     --cfg         Print scan configuration template"
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

case "$1" in
   --sch)
     if test -f $example_sched_file ; then
        cat $example_sched_file
        exit 0
     fi
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

# create temporary directory
out_dir="$(mktemp -d)"

# Get current time_t and convert to a UTC timestamp
timet=$(date +%s)
utc=$(date -u --date=@${timet} +%Y%m%dT%H%M%SZ)

plan_dirname="tempo_plan_${utc}"
target_dir="$out_dir/$plan_dirname"

if ! test -d $target_dir ; then
   mkdir -p $target_dir || error_exit "mkdir failed"
fi

echo "Writing temporary output to $target_dir"

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
   echo "WARNING: Scan tailoring file not specified"
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
  tar czf $tarfile --remove-files -C $(dirname $target_dir) ${plan_dirname}
  if test -d $out_dir ; then
    /bin/rmdir $out_dir
  fi
fi

echo "Created $tarfile"
