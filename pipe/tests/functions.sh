#
# Usage: wait_for_empty_dir $DIR $num_tries
#
# The wait time starts at 1 sec
# and doubles with each iteration.
#
wait_for_empty_dir()
{
  dir=$1
  n=$2
  t=1
  while test $n -ge 0; do
    sleep $t
    num_files=`ls -1 "$dir" | wc -l`
    if test $num_files -eq 0; then
       break
    fi
    n=`expr $n - 1`
    t=`expr $t \* 2`
  done
}
