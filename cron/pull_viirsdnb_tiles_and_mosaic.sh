#! /bin/sh

if test $# -lt 3 ; then
   echo "Usage: $(basename $0) destdir <source-url> <token-file> [date-args]"
   exit 0
fi

rootdir="$1"
source_url="$2"
token_file="$3"

# Remaining args specify which month, otherwise use the current month:
if test $# -gt 3 ; then
   shift 3
   yyyy_mm="$1"
else
   yyyy_mm="$(date --date today +%Y/%m)"
fi

if ! test -d "$rootdir" ; then
   mkdir -p $rootdir || exit 1
fi

incoming_dir="$rootdir/incoming"
mosaic_dir="$rootdir/mosaics"

# VIIRS-DNB monthly composites are labeled by the day-of-year number
# for the first day of each month:
yyyy_doy="$(date --date ${yyyy_mm}-01 +%Y/%j)"

# Exit silently if we already have a mosaic for this month:
tag="$(echo $yyyy_doy | tr -d /)"
mosaic_path="$(printf $mosaic_dir/VNP46A3_A%s.nc $tag)"
if test -f "$mosaic_path" ; then
   exit 0
fi

# Prepare to generate the mosaic:
if ! test -d $incoming_dir ; then
   mkdir -p $incoming_dir || exit 1
fi

if ! test -d $mosaic_dir ; then
   mkdir -p $mosaic_dir || exit 1
fi

#DRYRUN="--dryrun"
DRYRUN=""

log_message()
{
   _msg="$1"
   echo "$(date -u +%Y%m%d%H%M%SZ): $_msg"
}

# Download composite tiles to an 'incoming' directory:
log_message "checking for new VIIRS-DNB tiles: ${yyyy_doy}"
pull_viirsdnb_tiles.py $DRYRUN -t "$token_file" -s "${source_url}/${yyyy_doy}" -d $incoming_dir
if test "$?" -ne 0 ; then
   #log_message "*** Error:  VIIRS-DNB tile download failed: $yyyy_doy"
   exit 1
fi

num_incoming=$(find $incoming_dir -mindepth 1 -maxdepth 1 -type f | wc -l)
if test $num_incoming -eq 0 ; then
   log_message "no files downloaded"
   exit 0
fi

# Merge the tiles into a mosaic, storing the result in a separate directory:
log_message "generating VIIRS-DNB mosaic: ${yyyy_doy}"
srun --job-name DNB --nodes=1 --ntasks=1 \
     make_viirsdnb_mosaic.sh $incoming_dir $mosaic_dir
if test "$?" -ne 0 ; then
   log_message "*** Error:  VIIRS-DNB mosaic generation failed: $yyyy_doy"
   exit 1
fi

# Delete all tiles in the 'incoming' directory
find $incoming_dir -mindepth 1 -maxdepth 1 -type f -name "VNP46A3.*.h5" -delete
log_message "finished VIIRS-DNB mosaic: ${yyyy_doy}"

