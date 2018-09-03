#! /bin/sh

file_path="$1"

# move the file here for processing
/bin/mv "$file_path" .  || exit 1

# delete the rename directory which should now
# be empty
dir=`dirname "$file_path"`
/bin/rmdir "$dir" || exit 1
