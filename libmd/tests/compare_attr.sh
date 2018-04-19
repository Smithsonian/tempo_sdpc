#!/bin/sh

# Compare attributes in test_attr.nc with inputs
OTS_ROOT=`awk '{if($1=="OTS_ROOT")print $3}' < ../../Makefile.defines`

ncfile=test_attr.nc
nmlfile=boilerplate.nml

$OTS_ROOT/bin/ncdump -h $ncfile >| nc_header.txt

#comparison values
cmlat_in="0"
cmlon_in="0"
polylat_in="-10, -10, 10, 10"
polylon_in="-10, 10, 10, -10"
polyseq_in="1, 2, 3, 4"
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

exitstat=0
#compare
if [ "$cmlat_out" != "$cmlat_in" ] ; then exitstat=1 ; fi
if [ "$cmlon_out" != "$cmlon_in" ] ; then exitstat=2 ; fi
if [ "$polylat_out" != "$polylat_in" ] ; then exitstat=3 ; fi
if [ "$polylon_out" != "$polylon_in" ] ; then exitstat=4 ; fi
if [ "$polyseq_out" != "$polyseq_in" ] ; then exitstat=5 ; fi
if [ "$inputs_out" != "$inputs_in" ] ; then exitstat=6 ; fi
if [ "$local_gran_out" != "$local_gran_in" ] ; then exitstat=7 ; fi
if [ "$local_vers_out" != "$local_vers_in" ] ; then exitstat=8 ; fi
if [ "$colname_out" != "$colname_in" ] ; then exitstat=9 ; fi
if [ "$colvers_out" != "$colvers_in" ] ; then exitstat=10 ; fi
if [ "$platform_out" != "$platform_in" ] ; then exitstat=11 ; fi
if [ "$inst_out" != "$inst_in" ] ; then exitstat=12 ; fi
if [ "$acdesc_out" != "$acdesc_in" ] ; then exitstat=13 ; fi
if [ "$acval_out" != "$acval_in" ] ; then exitstat=14 ; fi
if [ "$abs_out" != "$abs_in" ] ; then exitstat=15 ; fi
if [ "$keys_out" != "$keys_in" ] ; then exitstat=16 ; fi

if [ $exitstat -ne 0 ]
then
  echo "*** test_attr: comparison failed"
fi

exit $exitstat
