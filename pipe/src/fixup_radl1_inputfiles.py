#! /usr/bin/env python3

import sys
import os
import re
import argparse

from netCDF4 import Dataset as NetCDFFile

def fix_input_files_entry (metfile, input_files):

    if not os.path.isfile (metfile):
        return

    with open (metfile, 'r') as fp:
        text = fp.read()

    text = text.replace ("${input_files}", "(" + ", ".join (f'"{w}"' for w in input_files) + ")")
    text = text.replace ("${#input_files}", "%d" % (len(input_files)))

    with open (metfile, 'w') as fp:
        fp.write(text)

def main ():
    parser = argparse.ArgumentParser(description='if necessary, expand input_files variable')
    parser.add_argument('radiance_file', help="Level 1 radiance file path")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    metfile = args.radiance_file + ".met"

    with NetCDFFile(args.radiance_file, "r") as nc:
        input_files = nc.getncattr ("input_files")

    fix_input_files_entry (metfile, input_files)

if __name__ == '__main__':
    main()
