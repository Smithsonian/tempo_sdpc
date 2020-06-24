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

def copy_frame (index, srcpath, destpath):
    src = Dataset (srcpath, 'r')

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

    uv_var = uv.variables[varname]
    dims = uv_var.shape
    if (dims[0] <= index):
        print ('copy_level1_frame.py: no frame with index={} ({} dims={})'.format(index, varname, dims))
        sys.exit(0)

    dest = Dataset (destpath, "w", format="NETCDF4")
    dest.setncattr ('source_file', srcpath)
    dest.setncattr ('index', index)

    slab = uv.variables[varname][index, :,:]
    define_target_group (dest, uv_group, varname, slab)

    vis = src.groups[vis_group]
    slab = vis.variables[varname][index, :,:]
    define_target_group (dest, vis_group, varname, slab)

    dest.close()
    src.close()

def main():
    parser = argparse.ArgumentParser(description='extract spectral arrays for one frame')
    parser.add_argument('index', help="frame index")
    parser.add_argument('infile', help="netCDF4 Level 1 file name")
    parser.add_argument('outfile', help="netCDF4 output file name")
    args = parser.parse_args()

    copy_frame (int(args.index), args.infile, args.outfile)

if __name__ == '__main__':
    main()
