#! /bin/sh

if test $# -le 1 ; then
  echo "Usage:  $0 ext <src-files>"
  exit 0
fi

ext=$1
shift
srcfiles=$@

mkfortranrule(){
  file=$1
  base=$(basename $file .${ext})

  mod_deps=""

  use_deps=`grep -o -i "^[ \t]*use[ \t]*[^ ,:]*" $file 2> /dev/null | sed s,use,,Ig | tr -d '\n'`
  if test "X${use_deps}" != "X" ; then
     for u in $use_deps; do
        mod_deps="\$(ODIR)/${u}.o $mod_deps"
     done
     mod_deps=`echo $mod_deps | tr ' ' '\n' | sort | uniq | grep -v ISO_C_BINDING | grep -v EZspline | grep -v HDF5 | tr '\n' ' '`
  fi

cat <<EOM
\$(ODIR)/${base}.o:	\$(SDIR)/${base}.${ext} ${mod_deps}
	\$(FC) -c \$(FCFLAGS) \$(IFLAGS) \$(SDIR)/${base}.${ext} -o \$(ODIR)/${base}.o \$(MODOUTFLG) \$(ODIR)
EOM
}

mkcrule(){
  file=$1
  base=$(basename $file .${ext})
cat <<EOM
\$(ODIR)/${base}.o:	\$(SDIR)/${base}.${ext}
	\$(CC) -c \$(CCFLAGS) \$(CINCFLGS) \$(SDIR)/${base}.${ext} -o \$(ODIR)/${base}.o
EOM
}

case $ext in
 f90 )
      for f in $srcfiles ; do
         mkfortranrule $f
      done
      ;;
   c )
      for f in $srcfiles ; do
         mkcrule $f
      done
      ;;
esac
