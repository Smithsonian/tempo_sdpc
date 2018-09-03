#! /bin/sh

: ${SDPC_ARCHIVE_DIR?}

granule_basename="$1"

case "${granule_basename}" in
  *irr0*)
     IFS='_' read -r prefix date version suffix<<<$granule_basename
     ;;

  *rad0*)
     IFS='_' read -r prefix date scan_num granule_num version suffix<<<$granule_basename
     ;;
esac

# Split date=YYYYMMDDTHHMMSSZ -> YYYY/MM/DD = date_subdirs
ymd=$(echo $date | cut -f 1 -d T)
date_subdirs=$(date --date $ymd +'%Y/%-m/%-d')

# Trim 'v' from version string
version_num=$(echo $version | tr -d v)

# tempo_*drk0.nc = unprocessed level 0 dark granule
# tempo_*drk1.nc = processed dark granule
dir_path="${SDPC_ARCHIVE_DIR}/L0/${version_num}/drk/${date_subdirs}"
if ! test -d "$dir_path" ; then
  echo $dir_path
else
  files=$(/bin/ls ${dir_path}/*drk1.nc)
  echo $files
fi

