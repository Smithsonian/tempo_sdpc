#! /bin/sh

file_path="$1"

# move the file here for processing
/bin/mv "$file_path" .  || exit 1

# do some nominal processing:
bname=`basename "$file_path"`
tar tf "$bname" | sort > tar.lst
chksum=`md5sum tar.lst | awk '{ print $1 }'`
touch "$chksum"

# delete the rename directory which should now
# be empty
dir=`dirname "$file_path"`
/bin/rmdir "$dir" || exit 1
