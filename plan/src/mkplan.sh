#! /bin/sh

PGMNAME="$(basename $0)"

if ! test $# -eq 1 ; then
   echo "Usage: $PGMNAME <schedule-config> | [options]"
   echo "   Options:"
   echo "     --sch     Print schedule configuration template"
   echo "     --cfg     Print scan configuration template"
   exit 0
fi

# SDPC_OTS_ROOT is used in the plan.cfg file
: "${SDPC_OTS_ROOT:?SDPC_OTS_ROOT is not set}"
: "${SDPC_ROOT:?SDPC_ROOT is not set}"

default_config_file="$SDPC_ROOT/share/plan.cfg"
example_sched_file="$SDPC_ROOT/share/plan_sched.sh.example"

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
   *)
     sched_file="$1"
     ;;
esac

if test -r "$sched_file" ; then
   . $sched_file
else
   echo "*** Error: cannot read schedule file: $sched_file"
   exit 1
fi

if test x"$prefix" = x ; then
   prefix="tempo_plan_$(date -u +%Y%m%d%H%M%SZ)"
fi

if test x"$config_file" = x ; then
   config_file="$default_config_file"
fi

: "${maneuver_file:?maneuver_file is not defined}"
: "${plan_type:?plan_type is not defined}"
: "${plan_num_days:?plan_num_days is not defined}"
: "${plan_start_day:?plan_start_day is not defined}"

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

echo "Writing temporary output to $target_dir"

_tailor="${target_dir}/${prefix}_scantailor.csv"
_master="${target_dir}/${prefix}_masterscan_tbl.csv"
_plan="${target_dir}/${prefix}_earthscan_tbl.csv"

# generate master scan table, and scan plan
plan -c $config_file -M $maneuver_file \
     -s $plan_type -d $plan_start_day -n $plan_num_days \
     $plan_options -o $_plan \
     -m $_master || error_exit "failed generating scan plan"

/bin/cp $config_file $target_dir/${prefix}_plan.cfg
/bin/cp $maneuver_file $target_dir

# generate scan tailoring file
if test x"$tailoring_file" = x ; then
: "${tailor_epoch:?tailor_epoch is not defined}"
: "${tailor_start_day:?tailor_start_day is not defined}"
: "${tailor_num_days:?tailor_num_days is not defined}"
   plan -c $config_file --epoch $tailor_epoch \
        -d $tailor_start_day -n $tailor_num_days \
        -T $_tailor || error_exit "failed generating scan tailoring file"
elif test -r "$tailoring_file" ; then
  /bin/cp $tailoring_file $_tailor
else
  error_exit "cannot read $tailoring_file"
fi

tarfile="${prefix}.tar.gz"

if test -d $target_dir ; then
  tar czf $tarfile --remove-files -C $(dirname $target_dir) ${prefix}
  if test -d $out_dir ; then
    /bin/rmdir $out_dir
  fi
fi

echo "Created $tarfile"
