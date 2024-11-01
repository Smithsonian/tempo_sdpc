#! /usr/bin/env python3

import sys
import numpy as np
import argparse
from netCDF4 import Dataset

def copy_geolocation_vars (g_src, g_dst):
    varlist_2d = ["latitude", "longitude",
                  "solar_azimuth_angle", "solar_zenith_angle",
                  "viewing_azimuth_angle", "viewing_zenith_angle"]
    varlist_3d = ["latitude_bounds", "longitude_bounds"]

    for var in varlist_2d:
        g_dst.variables[var][:,:] = g_src.variables[var][:,:]
    for var in varlist_3d:
        g_dst.variables[var][:,:,:] = g_src.variables[var][:,:,:]

def main():
    parser = argparse.ArgumentParser(description='Copy geolocation variables from RAD_L1 file to L2 file')
    parser.add_argument('radpath', help="Source netCDF4 RAD_L1 file name")
    parser.add_argument('l2path', help="Destination netCDF4 L2 file name")
    args = parser.parse_args()

    with Dataset (args.radpath, 'r') as src, Dataset (args.l2path, 'r+') as dst:
        copy_geolocation_vars (src.groups["band_290_490_nm"], dst.groups["geolocation"])

if __name__ == '__main__':
    main()
