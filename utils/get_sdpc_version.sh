#! /bin/sh

VN_DEF=4.7.0

if test -d .git ; then
   VN=$(git describe --tag  HEAD 2> /dev/null) || VN=$VN_DEF
else
   VN=$VN_DEF
fi

printf "$VN\n"
