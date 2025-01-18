#! /bin/sh

# configure test environment to get the right python installation
. ../etc/sdpc_env.sh

if test "$#" -ne 1 ; then
   echo "Usage: $(basename $0) register_py"
   exit 1
fi

register_py=$1

PROGNAME="$(basename $0)"
error_exit()
{
   echo "${PROGNAME}: ERROR: ${1:-'Unknown Error'}" 1>&2
   exit 1
}
trap error_exit ERR

testdir=$(mktemp -d)
dbfile="$testdir/iers.sqlite"

cleanup()
{
   if test -d "$testdir" ; then
      find $testdir -mindepth 1 -delete -name "bulletina*.txt"
      /bin/rm -f $dbfile
      /bin/rmdir $testdir
   fi
}

create_iers_file()
{
   line7="$1"
   filename=$(mktemp -u -p $testdir)

   cat <<- EOF > $filename
	line 0
	line 1
	line 2
	*IERSBULLETIN-A*
	line 4
	line 5
        line 6
	$line7
	EOF

   $register_py --dbfile $dbfile --rename $filename || error_exit "$register_py failed"
}

try_lookup()
{
   time="$1"
   expected="$2"
   if test "$#" -eq 3 ; then
       fail_expected="$3"
   else
       fail_expected="no"
   fi
   filename=$(mktemp -u -p $testdir)
   ./create_dummy_level1.py --tstart $time $filename

   got=$(../src/select_iers.py --dbfile $dbfile $filename)

   if test -f $filename ; then
      /bin/rm $filename
   fi

   if ! test x"$got" = x"$testdir/$expected" ; then
      if test x"$fail_expected" = xyes ; then
         echo "PASS: (lookup failed as expected)"
      else
         error_exit "FAIL: got:$got  expected:$testdir/$expected"
      fi
   else
      echo "PASS: $got"
   fi
}

echo "Testing IERS bulletin A file lookup:"

create_iers_file " 1 January 2013  Vol. x No. 1"
create_iers_file " 8 January 2013  Vol. x No. 1"
create_iers_file "15 January 2013  Vol. x No. 1"
create_iers_file " 1 March   2013  Vol. x No. 1"
create_iers_file " 6 June    2013  Vol. x No. 1"

# check day of update
try_lookup "2013-01-08T12:00:00Z" "bulletina_2013001_x_1.txt"

# check day before update
try_lookup "2013-01-07T12:00:00Z" "bulletina_2013001_x_1.txt"

# check day of update
try_lookup "2013-01-08T12:00:00Z" "bulletina_2013001_x_1.txt"

# check middle of gap
try_lookup "2013-02-28T12:00:00Z" "bulletina_2013015_x_1.txt"

# check day of update
try_lookup "2013-03-01T12:00:00Z" "bulletina_2013015_x_1.txt"
try_lookup "2013-06-10T12:00:00Z" "bulletina_2013157_x_1.txt"

# check failure
try_lookup "2013-12-15T12:00:00Z" "bulletina_2013157_x_1.txt" yes

cleanup

