#! /usr/bin/env python3

import os, sys
import glob
import argparse

from netCDF4 import Dataset as NetCDFFile

def classify_radiance_file (filename):
    with NetCDFFile (filename, "r") as nc:
        product_type = nc.product_type
        tstart = nc.time_coverage_start_since_epoch
        granule_flag = nc.variables["granule_flag"][0]
    info = {}
    info["tstart"] = tstart
    info["is_telem_only"] = product_type == 'RAD' and granule_flag == 4
    return info

def main():
    parser = argparse.ArgumentParser(description='Find/delete telemetry-only radiance files older than a specified file')
    parser.add_argument ('dir', help='directory path')
    parser.add_argument ('--delete', action='store_true', help='Delete only when this option is present')
    parser.add_argument ('--before', metavar='FILE', help='filename (will be deleted if telemetry-only)')
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    before_file = classify_radiance_file (args.before)
    tmax = before_file["tstart"]

    pattern = os.path.join (args.dir, "TEMPO_RAD_L1_V??_20??????T??????Z_S???G??.nc")
    files = glob.glob(pattern)
    for f in files:
        info = classify_radiance_file (f)
        if tmax < info["tstart"]:
            break
        if info["is_telem_only"]:
            if args.delete:
                os.remove(f)
            else:
                print(f)

if __name__ == "__main__":
    main()
