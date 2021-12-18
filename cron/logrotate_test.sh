#! /bin/sh

LOGDIR="var/log"
LOGFILES="asdc_geoscf.log asdc_goes.log cleanup.log geoscf.log goes.log ims.log"

LOGROTATE_ARGS="-f"

log_fast_fwd()
{
   num_weeks=$1

   sec_into_future=$((7*86400*$num_weeks))
   now=$(date +%s)

   for f in $LOGFILES ; do
      logfile="$LOGDIR/$f"
      echo "appending to $f"
      dd if=/dev/zero of=$logfile bs=1k count=10k oflag=append conv=notrunc > /dev/null 2>&1
      if test $sec_into_future -gt 0 ; then
        when=$(($now + $sec_into_future))
        tstamp=$(date --date @$when +%Y%m%d%H%M)
        touch -t $tstamp $logfile
      fi
   done
}

if test -d $LOGDIR ; then
  echo "Delete (or move) existing $LOGDIR directory tree before running this"
  exit 1
fi

mkdir -p $LOGDIR || exit 1

weeks=$(seq 0 6)
this_dir=$(pwd)

for num_weeks in $weeks ; do
   log_fast_fwd $num_weeks
   echo "==== Before rotate (num_weeks=$num_weeks):"
   ls -l $LOGDIR
   /usr/sbin/logrotate $LOGROTATE_ARGS -s var/log/logstate $this_dir/logrotate.conf
   echo "==== After rotate (num_weeks=$num_weeks):"
   ls -l $LOGDIR
done
