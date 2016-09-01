#! /bin/sh

# Comment out include lines in specified fortran file.
# Why do we want to do this?  See doxy_filter_inc.sh
# for details.

sed "s,include ',\! include ',g" $1

