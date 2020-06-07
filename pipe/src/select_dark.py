#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import sqlite3
from datetime import date
from netCDF4 import Dataset as NetCDFFile
import numpy as np

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
    exposure_time = nc.variables["exposure_time"][:]
    num_coadds = nc.variables["num_coadds"][:]
    keys["mean_exposure_time_per_coadd"] = np.mean(exposure_time/num_coadds)
    nc.close()
    return keys

def select_matching_dark (c, keys):
    tstart = keys["time_coverage_start_since_epoch"]
    exposure_time_per_coadd = keys["mean_exposure_time_per_coadd"]
    # FIXME - we'll also need to filter this on start time, but that's tricky when the times in the
    # test files are incorrect...
    cmd = "select path from \
    (      select path, min(mean_exposure_time_per_coadd) as dt from 'DRK_L1' where mean_exposure_time_per_coadd >= {} \
    union select path, max(mean_exposure_time_per_coadd) as dt from 'DRK_L1' where mean_exposure_time_per_coadd < {}) \
    order by abs({}-dt);".format (exposure_time_per_coadd, exposure_time_per_coadd, exposure_time_per_coadd)
    c.execute(cmd)
    rows = c.fetchall()
    rows = [r for r in rows if None not in r]
    return rows

def main():

    if len (sys.argv) == 1:
        print('Usage:  select_dark L0_FILE')
        sys.exit(0)

    level0_file = sys.argv[1]
    level0_keys = get_level0_keys (level0_file)

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
        darks = select_matching_dark (c, level0_keys)

    best_match = darks[0]
    print(best_match[0])

if __name__ == "__main__":
    main()
