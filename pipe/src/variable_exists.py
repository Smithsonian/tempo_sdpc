#! /usr/bin/env python3

import os
import sys
import argparse

from netCDF4 import Dataset as NetCDFFile

def main():
    parser = argparse.ArgumentParser(description='determine if netcdf4 variable exists')
    parser.add_argument('--var', help="variable path")
    parser.add_argument('ncfile', help="netCDF data file name")
    args = parser.parse_args()

    path = args.var
    var = os.path.basename (path)

    with NetCDFFile (args.ncfile, "r") as nc:
        if '/' in path:
            grp = nc[os.path.dirname (path)]
        else:
            grp = nc
        if var in grp.variables:
            print("yes")
        else:
            print("no")

if __name__ == '__main__':
    main()
