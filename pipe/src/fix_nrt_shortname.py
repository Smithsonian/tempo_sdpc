#! /usr/bin/env python3

import sys
import os
import re
import argparse
from netCDF4 import Dataset as NetCDFFile

def process_file (ncfile):

    # Extract the shortname from the ncfile name
    basename = os.path.basename (ncfile)
    if '_NRT_' not in basename:
        return

    m = re.match (r'(TEMPO_\w+_L\d_NRT)_V', basename)
    shortname = m.group(1)

    # Set shortname in the ncfile header
    with NetCDFFile (ncfile, "r+") as nc:
        nc.setncattr ('shortname', shortname)

    # Set shortname in the .met file
    metfile = ncfile + ".met"

    with open (metfile, 'r') as fp:
        text = fp.read()

    def shortname_endswith_nrt (matchobj):
        s = matchobj.group(0)
        return s.replace (matchobj.group(1), shortname)

    text = re.sub (r'SHORTNAME\s+NUM_VAL\s+=\s+1\s+VALUE\s+=\s+"([^"]+)"', shortname_endswith_nrt, text)

    with open (metfile, 'w') as fp:
        fp.write(text)

def main ():
    parser = argparse.ArgumentParser(description='Ensure shortname attribute has _NRT as needed')
    parser.add_argument('ncfile', help="Path to netcdf4 data product file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if "_NRT_" in args.ncfile:
        process_file (args.ncfile)

if __name__ == '__main__':
    main()
