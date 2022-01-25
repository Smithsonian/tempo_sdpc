#! /bin/sh

export SDPC_REFDATA_DIR=/tempo/nas0/sdpc_soft/refdata
export SDPC_ANCILLARY_ROOT=/tempo/nas0/sdpc/ancillary

export SDPC_GEOSCF_CONFIG="clim_config.ini"

#PREFIX="gdb --args"
#PREFIX="valgrind --tool=memcheck --log-file=valgrind.log --leak-check=yes --error-limit=no --num-callers=25"
#PREFIX="valgrind --tool=memcheck --log-file=valgrind.log --leak-check=full --show-leak-kinds=all --error-limit=no --num-callers=25"
#PREFIX="/usr/bin/time -v"
#PREFIX=""

$PREFIX ./test_clim > out.dat
#$PREFIX ./test_clim

##
## Use h5dump output to check values from the point sampled by test_clim.f90
##
fcast_file="$SDPC_ANCILLARY_ROOT/var/geoscf/2022/025/GEOS-CF.v01.fcst.sat_inst_1hr_r720x361_v72.20220124_12z_20220125_1800z_reorder.nc4"
h5dump -w 1 -d T -A 0 -s 0,380,144,0 -c 1,1,1,72 $fcast_file > T_hr18_lon-85,lat+36.dat
h5dump -w 1 -d U2M -A 0 -s 0,380,144 -c 1,1,1 $fcast_file > U2M_hr18_lon-85,lat+36.dat
h5dump -w 1 -d V2M -A 0 -s 0,380,144 -c 1,1,1 $fcast_file > V2M_hr18_lon-85,lat+36.dat
