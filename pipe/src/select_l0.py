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
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    return db_file_path

def get_keys (filename):
    nc = NetCDFFile (filename, "r")
    keys = {}
    keys["time_coverage_start_since_epoch"] = nc.getncattr ('time_coverage_start_since_epoch')
    keys["time_coverage_end_since_epoch"] = nc.getncattr ('time_coverage_end_since_epoch')
    nc.close()
    return keys

def select_matching (c, table_name, tstart, tend):
    cmd = "select path from '{}' where {} < time_coverage_end_since_epoch AND \
                                    time_coverage_start_since_epoch < {};".format(table_name, tstart, tend)
    c.execute(cmd)
    rows = c.fetchall()
    rows = [r for r in rows if None not in r]
    return rows

def main():
    parser = argparse.ArgumentParser(description='Select the appropriate Level 0 files')
    parser.add_argument ('--wait', type=float, default=60.0, help="Optional wait time")
    parser.add_argument ('--table', help="Product table name")
    parser.add_argument ('--begin', type=float, help="start time")
    parser.add_argument ('--end', type=float, help="end time")
    parser.add_argument ('--window', metavar='SECONDS', default=300.0, type=float,
                         help="Time coverage padding in minutes")
    parser.add_argument ('--granule', help="Path to relevant product file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if (args.window <= 0):
        print ('*** {}: invalid window = {}'.format(os.path.basename(sys.argv[0]), args.window))
        sys.exit(1)

    window_sec = args.window

    if args.granule:
        keys = get_keys (args.granule)
        tbeg = keys["time_coverage_start_since_epoch"] - window_sec
        tend = keys["time_coverage_end_since_epoch"] + window_sec
    else:
        tbeg = args.begin
        tend = args.end

    if args.wait > 0:
        num_tries_remaining = 2
    else:
        num_tries_remaining = 1

    db_path = get_db_path()

    """
     For back-compatibility sqlite has foreign keys turned off by default,
     and foreign_keys=off is ALWAYS stored in the database, regardless of
     the runtime setting when the database was created.  For this reason,
     we apparently need to turn it on explicitly, each time the database
     connection is established.
    """

    while num_tries_remaining > 0:
        num_tries_remaining = num_tries_remaining - 1
        try:
            with sqlite3.connect (db_path) as conn:
                conn.execute("pragma foreign_keys=on")
                c = conn.cursor()
                files = select_matching (c, args.table, tbeg, tend)
        except:
            files = []
        if len(files) > 0:
            break
        if args.wait > 0:
            time.sleep (args.wait)

    if len(files) == 0:
        print("NONE")
    else:
        path_list = list(map (' '.join, files))
        basenames = [os.path.basename(p) for p in path_list]
        indices = np.argsort(basenames)
        for i in indices:
            print(path_list[i])

if __name__ == "__main__":
    main()
