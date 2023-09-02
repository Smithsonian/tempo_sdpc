#! /bin/sh

if test "$#" -lt 1 ; then
   echo "Usage:  $0 FILE [FILE] ..."
   exit 0
fi

outdir=deflated
if ! test -d $outdir ; then
    mkdir -p $outdir
fi

deflate_file()
{
  infile=$1

  outfile="$outdir/$(basename $infile)"

  # ABI band 1 center = 470 nm  => /band_290_490_nm/spectral_channel=901
  # ABI band 2 center = 640 nm  => /band_540_740_nm/spectral_channel=520

  if test -f $infile ; then
     echo "deflating $infile"
     ncks -d spectral_channel,901 -g band_290_490_nm $infile $outfile
     ncks -A -d spectral_channel,520 -g band_540_740_nm $infile $outfile
     ncks -A -g geometry $infile $outfile
     ncks -A -g inr_input $infile $outfile
     ncks -A -v granule_flag,exposure_time,earth_sun_distance $infile $outfile
     gzip $outfile
  else
     echo "*** File not found: $infile"
  fi
}

file_list="$@"

for file in $file_list ; do
   deflate_file $file
done

