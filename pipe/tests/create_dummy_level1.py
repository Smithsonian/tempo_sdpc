#! /usr/bin/env python3

import os, sys
import argparse
from netCDF4 import Dataset as NetCDFFile
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description='create dummy Level 1 file')
    parser.add_argument('ncfile', help="netCDF data file name", default=None)
    parser.add_argument('--tstart', help="UTC start time, YYYY-MM-DDTHH:MM:SSZ", default=None)
    if len(sys.argv) < 3:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    with NetCDFFile(args.ncfile, "w", format="NETCDF4") as nc:
        nc.setncattr ("time_coverage_start", args.tstart)

if __name__ == "__main__":
    main()

