#! /usr/bin/env python3

import sys
import argparse

from netCDF4 import Dataset as NetCDFFile

def main():
    parser = argparse.ArgumentParser(description='print radiance header attribute')
    parser.add_argument('--attr', help="attribute name")
    parser.add_argument('ncfile', help="netCDF data file name")
    args = parser.parse_args()

    with NetCDFFile (args.ncfile, "r") as nc:
        if args.attr in nc.__dict__:
            value = nc.__dict__[args.attr]
        else:
            value = None

    if value != None:
        print(value)
    else:
        print ("*** Nonexistent attribute: {}".format(args.attr))
        sys.exit(1)

if __name__ == '__main__':
    main()
