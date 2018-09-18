#! /bin/sh

: ${SDPC_ARCHIVE_DIR?}

granule_basename="$1"

# prefix=TEMPO
# ptype=irr | rad
# level=Ld
# version=Vdd
# date=YYYYMMDDTHHMMSSZ

case "$granule_basename" in
 *irr*)
     IFS='_' read -r prefix ptype level version date<<<$granule_basename
     ;;
 *rad*)
     #suffix = "SdddGdd"
     IFS='_' read -r prefix ptype level version date suffix<<<$granule_basename
     ;;
  *)
     echo "$0: Unexpected filename format: $granule_basename"
     exit 1
    ;;
esac

# Split date=YYYYMMDDTHHMMSSZ -> YYYY/MM/DD = date_subdirs
ymd=$(echo $date | cut -f 1 -d T)
date_subdirs=$(date --date $ymd +'%Y/%-m/%-d')

# Trim 'V' and leading zeros from version number
version_num_with_leading_zeros=$(echo $version | tr -d V)
version_num=$((10#${version_num_with_leading_zeros}))

dir_path="${SDPC_ARCHIVE_DIR}/L0/${version_num}/drk/${date_subdirs}"
if ! test -d "$dir_path" ; then
  echo $dir_path
else
  files=$(/bin/ls ${dir_path}/TEMPO_drkt_*.nc)
  echo $files
fi

