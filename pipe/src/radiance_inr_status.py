#! /usr/bin/env python3

import argparse

from netCDF4 import Dataset as NetCDFFile

def main():
    parser = argparse.ArgumentParser(description='check radiance file INR status')
    parser.add_argument('ncfile', help="netCDF data file name")
    args = parser.parse_args()

    with NetCDFFile (args.ncfile, "r") as nc:
        print (nc.inr_status)

if __name__ == '__main__':
    main()
