#! /usr/bin/env python3

import os
import sys
import argparse
import numpy as np

from netCDF4 import Dataset as NetCDFFile

def main():
    parser = argparse.ArgumentParser(description='check if exposure time per frame exceeds threshold')
    parser.add_argument('--threshold', default=0.2, type=float,
                         help="Exposure time per frame [sec]")
    parser.add_argument('--verbose', action='store_true',
                         help="Verbose")
    parser.add_argument('ncfile', help="netCDF data file name")
    args = parser.parse_args()

    with NetCDFFile (args.ncfile, "r") as nc:
        exptime = nc.variables['exposure_time'][:]
        num_coadds = nc.variables['num_coadds'][:]
        mean_exptime_per_frame = np.mean(exptime/num_coadds)

    if args.verbose:
        print('mean_exptime_per_frame = {}'.format(mean_exptime_per_frame))

    if mean_exptime_per_frame > args.threshold:
        result=1
    else:
        result=0

    print(result)

if __name__ == '__main__':
    main()
