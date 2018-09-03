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
   IFS='_' read -r prefix date scan_num granule_num version suffix<<<$granule_basename

   # Split date=YYYYMMDDTHHMMSSZ -> YYYY/MM/DD = date_subdirs
   ymd=$(echo $date | cut -f 1 -d T)
   date_subdirs=$(date --date $ymd +'%Y/%-m/%-d')

   # Trim 'v' from version string
   version_num=$(echo $version | tr -d v)

   # tempo_*irr1.nc = Level 1 irradiance
   dir_path="${SDPC_ARCHIVE_DIR}/L1/${version_num}/irr/${date_subdirs}"
   if ! test -d "$dir_path" ; then
      use_fixed_irr
   else
      files=$(/bin/ls ${dir_path}/*irr1.nc)
      if test x"$files" != x ; then
         echo $files
      else
         use_fixed_irr
      fi
   fi
}

lookup_irr "$granule_basename"
