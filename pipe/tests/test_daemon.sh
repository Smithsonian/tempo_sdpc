#! /bin/sh

/bin/rm -f hello.log

export SLANG_MODULE_PATH="../src"
./hello.sl || exit 1

sleep 2

num_hellos=`grep -c Hello hello.log`
if ! test "$num_hellos" -eq 3 ; then
   exit 1
fi

/bin/rm hello.log || exit 1
