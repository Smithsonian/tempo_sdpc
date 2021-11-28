#! /bin/sh

PGMNAME="mkplan.sh"

if ! test $# -eq 1 ; then
   echo "Usage: $PGMNAME <schedule-config>"
   exit 0
fi

sched_file="$1"

if test -r "$sched_file" ; then
   . $sched_file
else
   echo "*** Error: cannot read schedule file: $sched_file"
   exit 1
fi

: "${config_file:?config_file is not defined}"
: "${maneuver_file:?maneuver_file is not defined}"

: "${plan_type:?plan_type is not defined}"
: "${plan_num_days:?plan_num_days is not defined}"
: "${plan_start_day:?plan_start_day is not defined}"

: "${prefix:?prefix is not defined}"

# SDPC_OTS_ROOT is used in the plan.cfg file
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT is not set}"
: "${SDPC_ROOT:?SDPC_ROOT is not set}"

error_exit()
{
  echo "$PGMNAME [$$]: ERROR: ${1:-'Unknown Error'}" 1>&2
  exit 1
}
trap error_exit ERR

out_dir="$(mktemp -d)"
target_dir="$out_dir/$prefix"

if ! test -d $target_dir ; then
   mkdir -p $target_dir || error_exit "mkdir failed"
fi

echo "Writing output to $target_dir"

_tailor="${target_dir}/${prefix}_scantailor.csv"
_master="${target_dir}/${prefix}_masterscan_tbl.csv"
_plan="${target_dir}/${prefix}_earthscan_tbl.csv"

# generate master scan table, and scan plan
plan -v -c $config_file -M $maneuver_file \
     -s $plan_type -d $plan_start_day -n $plan_num_days \
     -o $_plan \
     -m $_master || error_exit "failed generating scan plan"

/bin/cp $config_file $target_dir/${prefix}_plan.cfg
/bin/cp $maneuver_file $target_dir

# generate scan tailoring file
if test x"$tailoring_file" = x ; then
   plan -c $config_file --epoch $tailor_epoch \
        -d $tailor_start_day -n $tailor_num_days \
        -T $_tailor || error_exit "failed generating scan tailoring file"
elif test -r "$tailoring_file" ; then
  /bin/cp $tailoring_file $_tailor
else
  error_exit "cannot read $tailoring_file"
fi

if test -d $target_dir ; then
  tar cvzf ${prefix}.tar.gz --remove-files -C $(dirname $target_dir) ${prefix}
  if test -d $out_dir ; then
    /bin/rmdir $out_dir
  fi
fi
