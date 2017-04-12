#! /bin/sh

ROOT_DIR="/vex/d2/tempo/plans/ancillary_data/modis_land_cover"

# MRT_DATA_DIR is required by 'mrt'
export MRT_DATA_DIR=${ROOT_DIR}/mrt/data

OUT_DIR=modis_files

VLIST="2 3 4 5 6 7"

create_mosaic(){

 __SPECTRAL_SUBSET="$1"

 for VV in $VLIST ; do
   files=$(ls modis_files/hdf/M*v0${VV}*.hdf)
   printf "%s\n" $files > modis_hdf.lis
   mrt/bin/mrtmosaic -i modis_hdf.lis -s "${__SPECTRAL_SUBSET}" \
         -o ${OUT_DIR}/hdf_mosaic/modis_v${VV}_mosaic.hdf
 done
}

reformat_hdf_to_tif(){

 __SPECTRAL_SUBSET="$1"
 SUBDIR=$2

# output projection parameter is radius of authalic sphere [meters]

cat << EOF > mrt_param.prm
INPUT_FILENAME = cmdline_override.hdf
OUTPUT_FILENAME = cmdline_override.tif
SPECTRAL_SUBSET = (${__SPECTRAL_SUBSET})
OUTPUT_PROJECTION_TYPE = GEO
OUTPUT_PROJECTION_PARAMETERS = (6371007.181000)
# OUTPUT_PIXEL_SIZE = 0.01
EOF

TARGET_DIR=$OUT_DIR/${SUBDIR}

for VV in $VLIST ; do
 infile=${OUT_DIR}/hdf_mosaic/modis_v${VV}_mosaic.hdf
 bn=$(basename $infile .hdf)
 mrt/bin/resample -p mrt_param.prm -i $infile \
      -o ${TARGET_DIR}/${bn}.tif
done
}

# select Land_Cover_Type_1 for land cover bits and
#    and Land_Cover_Type_QC for land/water mask
SPECTRAL_SUBSET="1 0 0 0 0 0 0 0 0 0 1 0 0 0 0"

# 1. create latitude band mosaics in HDF format
create_mosaic "$SPECTRAL_SUBSET"

# 2. resample HDF to geotiff:
reformat_hdf_to_tif "1 1" tif_mosaic
