#! /usr/bin/env python3

import sys
import numpy as np
import argparse
from netCDF4 import Dataset

def define_target_group (dest, grpname, varname, cube):
    dims = cube.shape
    grp = dest.createGroup (grpname)
    grp.createDimension ('xtrack', dims[0])
    grp.createDimension ('spectral_channel', dims[1])
    var_rad = grp.createVariable (varname, 'f4', ('xtrack', 'spectral_channel'))
    var_rad[:,:] = cube[:,:]

def copy_mirror_step (step, srcpath, destpath):
    src = Dataset (srcpath, 'r')
    mirror_step = src.variables['mirror_step'][:]
    result = np.where (mirror_step == step)
    if len(result[0]) == 0:
        print('*** Error: mirror step {} not found, file = {}'.format(step, srcpath))
        sys.exit(1)
    k = result[0][0]

    uv_group = 'band_290_490_nm'
    vis_group = 'band_540_740_nm'

    uv = src.groups[uv_group]

    if 'radiance' in uv.variables.keys():
        varname = 'radiance'
    elif 'irradiance' in uv.variables.keys():
        varname = 'irradiance'
    else:
        print ('*** Error: unsupported file: {}'.format(srcpath))
        sys.exit(1)

    dest = Dataset (destpath, "w", format="NETCDF4")
    dest.setncattr ('source_file', srcpath)
    dest.setncattr ('mirror_step', step)
    slab = uv.variables[varname][k, :,:]
    define_target_group (dest, uv_group, varname, slab)

    vis = src.groups[vis_group]
    slab = vis.variables[varname][k, :,:]
    define_target_group (dest, vis_group, varname, slab)

    dest.close()
    src.close()

def main():
    parser = argparse.ArgumentParser(description='extract spectral arrays for one mirror step')
    parser.add_argument('mirror_step', help="mirror_step index")
    parser.add_argument('infile', help="netCDF4 Level 1 file name")
    parser.add_argument('outfile', help="netCDF4 output file name")
    args = parser.parse_args()

    copy_mirror_step (int(args.mirror_step), args.infile, args.outfile)

if __name__ == '__main__':
    main()
