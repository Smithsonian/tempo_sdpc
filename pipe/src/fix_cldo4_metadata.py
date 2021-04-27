#! /usr/bin/env python3

import os, sys
import argparse
from netCDF4 import Dataset as NetCDFFile

def main():
    parser = argparse.ArgumentParser(description='fix CLDO4 product metadata')
    parser.add_argument('ncfile', help="netCDF data file name", default=None)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    with NetCDFFile(args.ncfile, "r+") as nc:
        nc.setncattr ("product_type", "CLDO4")

if __name__ == "__main__":
    main()

