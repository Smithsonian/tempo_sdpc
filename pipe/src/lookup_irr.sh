#! /bin/sh

: ${SDPC_RUN_DIR?}

granule_basename="$1"

use_fixed_irr()
{
  irr_basename=`cat ${SDPC_RUN_DIR}/L1/irradiance/out/IRR_LOOKUP_PLACEHOLDER`
  echo "${SDPC_RUN_DIR}/L1/irradiance/out/${irr_basename}.nc"
}

lookup_irr()
{
   IFS='_' read -r prefix ptype level version date suffix<<<$granule_basename

   # Split date=YYYYMMDDTHHMMSSZ -> YYYY/MM/DD = date_subdirs
   ymd=$(echo $date | cut -f 1 -d T)
   date_subdirs=$(date --date $ymd +'%Y/%-m/%-d')

   # Trim 'V' and leading zeros from version number
   version_num_with_leading_zeros=$(echo $version | tr -d V)
   version_num=$((10#${version_num_with_leading_zeros}))

   dir_path="${SDPC_ARCHIVE_DIR}/L1/${version_num}/irr/${date_subdirs}"
   if ! test -d "$dir_path" ; then
      use_fixed_irr
   else
      files=$(/bin/ls ${dir_path}/TEMPO_irr_*.nc)
      if test x"$files" != x ; then
         echo $files
      else
         use_fixed_irr
      fi
   fi
}

lookup_irr "$granule_basename"
