#! /usr/bin/env python3

# for eprint definition
#from __future__ import print_function

import os, sys
import time
import sqlite3
from datetime import date
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import argparse

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_db_path():
    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if (arch_dir == None):
        eprint ('*** Error: SDPC_ARCHIVE_DIR is not set')
        sys.exit(1)

    db_basename = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if (db_basename == None):
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    db_dir = os.path.join (arch_dir, "registry")
    db_path = os.path.join (db_dir, db_basename)

    return db_path

def get_level0_keys (filename):
    nc = NetCDFFile (filename, "r")
    keys = {}
    keys["time_coverage_start_since_epoch"] = nc.getncattr ('time_coverage_start_since_epoch')
    keys["time_coverage_end_since_epoch"] = nc.getncattr ('time_coverage_end_since_epoch')
    nc.close()
    return keys

def select_matching_hk (c, window_min, keys):
    tpad = window_min * 60.0
    tstart = keys["time_coverage_start_since_epoch"] - tpad
    tend = keys["time_coverage_end_since_epoch"] + tpad
    cmd = "select path from 'HK_L0' where {} < time_coverage_end_since_epoch AND \
                                    time_coverage_start_since_epoch < {};".format(tstart, tend)
    c.execute(cmd)
    rows = c.fetchall()
    rows = [r for r in rows if None not in r]
    return rows

def main():
    parser = argparse.ArgumentParser(description='Select the appropriate Level 0 HK files')
    parser.add_argument ('--window', metavar='MINUTES', default=2.0, type=float,
                         help="Time coverage padding in minutes")
    parser.add_argument ('level0_file', help="Path to relevant Level 0 file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if (args.window <= 0):
        print ('*** {}: invalid window = {}'.format(os.path.basename(sys.argv[0]), args.window))
        sys.exit(1)

    level0_keys = get_level0_keys (args.level0_file)

    db_path = get_db_path()

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")
    with conn:
        c = conn.cursor()
        hks = select_matching_hk (c, args.window, level0_keys)

    if len(hks) == 0:
        print("NONE")
    else:
        path_list = list(map (' '.join, hks))
        basenames = [os.path.basename(p) for p in path_list]
        indices = np.argsort(basenames)
        for i in indices:
            print(path_list[i])

if __name__ == "__main__":
    main()
