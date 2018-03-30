#!/bin/sh

# Compare attributes in test_attr.nc with inputs
OTS_ROOT=`awk '{if($1=="OTS_ROOT")print $3}' < ../../Makefile.defines`

ncfile=test_attr.nc
nmlfile=boilerplate.nml

$OTS_ROOT/bin/ncdump -h $ncfile >| nc_header.txt

#comparison values
cmlat_in="-1"
cmlon_in="-2"
polylat_in="-10, -10, -10, -10, -10, -5, 0, 5, 10, 10, 10, 10, 10, 5, 0, -5"
polylon_in="-10, -6, 0, 4, 10, 10, 10, 10, 10, 4, 0, -6, -10, -10, -10, -10"
polyseq_in="1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16"
inputs_in="l1_radiance.nc, l1_irradiance.nc, lookup.txt, reference_table.nc, otherstuff.dat"
local_gran_in=$ncfile
local_vers_in="(1)"
colname_in=`awk -F\' '/collection_shortname/{print $2}' boilerplate.nml`
colvers_in=`awk -F\' '/collection_version/{print $2}' boilerplate.nml`
platform_in=`awk -F\' '/platform/{print $2}' boilerplate.nml`
inst_in="TEMPO"
acdesc_in=`awk -F\' '/access_description/{print $2}' boilerplate.nml`
acval_in=`awk -F\' '/access_value/{print $2}' boilerplate.nml`
abs_in=`awk -F\" '/abstract/{print $2}' boilerplate.nml`
keys_in=`awk -F\" '/keywords/{print $2}' boilerplate.nml`

#attribute values
cmlat_out=`cat nc_header.txt | awk -F=\  '/:centroid_mean_latitude/{print $2}' | awk -F. '{print $1}'`
cmlon_out=`cat nc_header.txt | awk -F=\  '/:centroid_mean_longitude/{print $2}' | awk -F. '{print $1}'`
polylat_out=`cat nc_header.txt | awk -F=\  '/:polygon_latitude/{print $2}' | sed s/.f//g | sed s/\ \;//`
polylon_out=`cat nc_header.txt | awk -F=\  '/:polygon_longitude/{print $2}' | sed s/.f//g | sed s/\ \;//`
polyseq_out=`cat nc_header.txt | awk -F=\  '/:polygon_sequence/{print $2}' | sed s/\ \;//`
inputs_out=`cat nc_header.txt | awk -F\" '/:input_files/{print $2}'`
local_gran_out=`cat nc_header.txt | awk -F\" '/:local_granule_id/{print $2}'`
local_vers_out=`cat nc_header.txt | awk -F\" '/:local_version_id/{print $2}'`
colname_out=`cat nc_header.txt | awk -F\" '/:collection_shortname/{print $2}'`
colvers_out=`cat nc_header.txt | awk -F\" '/:collection_version/{print $2}'`
platform_out=`cat nc_header.txt | awk -F\" '/:platform/{print $2}'`
inst_out=`cat nc_header.txt | awk -F\" '/:instrument/{print $2}'`
acdesc_out=`cat nc_header.txt | awk -F\" '/:access_description/{print $2}'`
acval_out=`cat nc_header.txt | awk -F\" '/:access_value/{print $2}'`
abs_out=`cat nc_header.txt | awk -F\" '/:abstract/{print $2}' | sed "s%[\\]%%g"`
keys_out=`cat nc_header.txt | awk -F\" '/:keywords/{print $2}'`

problem=0
#compare
if [ "$cmlat_out" != "$cmlat_in" ] ; then problem=1 ; fi
if [ "$cmlon_out" != "$cmlon_in" ] ; then problem=1 ; fi
if [ "$polylat_out" != "$polylat_in" ] ; then problem=1 ; fi
if [ "$polylon_out" != "$polylon_in" ] ; then problem=1 ; fi
if [ "$polyseq_out" != "$polyseq_in" ] ; then problem=1 ; fi
if [ "$inputs_out" != "$inputs_in" ] ; then problem=1 ; fi
if [ "$local_gran_out" != "$local_gran_in" ] ; then problem=1 ; fi
if [ "$local_vers_out" != "$local_vers_in" ] ; then problem=1 ; fi
if [ "$colname_out" != "$colname_in" ] ; then problem=1 ; fi
if [ "$colvers_out" != "$colvers_in" ] ; then problem=1 ; fi
if [ "$platform_out" != "$platform_in" ] ; then problem=1 ; fi
if [ "$inst_out" != "$inst_in" ] ; then problem=1 ; fi
if [ "$acdesc_out" != "$acdesc_in" ] ; then problem=1 ; fi
if [ "$acval_out" != "$acval_in" ] ; then problem=1 ; fi
if [ "$abs_out" != "$abs_in" ] ; then problem=1 ; fi
if [ "$keys_out" != "$keys_in" ] ; then problem=1 ; fi

if [ $problem -ne 0 ]
then
  echo "*** test_attr: comparison failed"
  exitstat=1
fi

exit $exitstat
