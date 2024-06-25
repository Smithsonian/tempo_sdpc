#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

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

def get_level0_keys (filename):
    nc = NetCDFFile (filename, "r")
    keys = {}
    keys["time_coverage_start_since_epoch"] = nc.getncattr ('time_coverage_start_since_epoch')
    exposure_time = nc.variables["exposure_time"][:]
    num_coadds = nc.variables["num_coadds"][:]
    keys["mean_exposure_time_per_coadd"] = np.mean(exposure_time/num_coadds)
    nc.close()
    return keys

def select_matching_dark (c, window_hours, keys):
    dt_name="mean_exposure_time_per_coadd"
    t_name="time_coverage_start_since_epoch"

    tx = keys[t_name]
    dtx = keys[dt_name]
    t_ok = window_hours * 3600.0   # time offset acceptable for dark measurements

    subst = {'dtx':dtx, 'tx':tx, 'beg':tx-t_ok, 'end':tx+t_ok}

    cmd = "select path from DRK_L1 where istart between :beg and :end order by abs({dt_name}-:dtx),abs(istart-:tx) limit 1".format(**locals())
    c.execute(cmd, subst)

    rows = c.fetchone()
    if len(rows) > 0:
        return rows[0]
    else:
        return ""

def main():
    parser = argparse.ArgumentParser(description='Select the appropriate Level 1 dark file')
    parser.add_argument ('--window', metavar='HOURS', default=36.0, type=float,
                         help="Acceptable time offset in hours")
    parser.add_argument ('level0_file', help="Path to Level 0 file needing dark correction")
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

    with sqlite3.connect ("file:{}?mode=ro".format(db_path), uri=True, timeout=20.0) as conn:
        conn.execute("pragma foreign_keys=on")
        #conn.set_trace_callback(print)
        c = conn.cursor()
        path = select_matching_dark (c, args.window, level0_keys)

    print(path)

if __name__ == "__main__":
    main()
