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

def get_file_keys (filename):
    nc = NetCDFFile (filename, "r")
    keys = {}
    keys["time_coverage_start_since_epoch"] = nc.getncattr ('time_coverage_start_since_epoch')
    nc.close()
    return keys

def select_irradiance (c, window_days, keys):
    t_name="time_coverage_start_since_epoch"

    tx = keys[t_name]
    t_ok = window_days * 86400.0;  # max time offset for irradiance measurement

    subst = {'tx':tx, 'beg':tx-t_ok, 'end':tx}

    cmd = "select path from 'IRR_L1' \
         where ({t_name} > :beg and {t_name} < :end) order by abs({t_name}-:tx);".format (**locals())

    c.execute(cmd, subst)
    rows = c.fetchall()
    rows = [r for r in rows if None not in r]

    return rows

def main():
    parser = argparse.ArgumentParser(description='Select the appropriate Level 1 irradiance file')
    parser.add_argument ('--window', metavar='DAYS', default=7.0, type=float,
                         help="Acceptable time offset in days")
    parser.add_argument ('level1_file', help="Path to Level 1 radiance file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if (args.window <= 0):
        print ('*** {}: invalid window = {}'.format(os.path.basename(sys.argv[0]), args.window))
        sys.exit(1)

    file_keys = get_file_keys (args.level1_file)

    db_path = get_db_path()

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    with sqlite3.connect (db_path) as conn:
        conn.execute("pragma foreign_keys=on")
        #conn.set_trace_callback(print)
        c = conn.cursor()
        irr = select_irradiance (c, args.window, file_keys)

    best_match = irr[0]
    print(best_match[0])

if __name__ == "__main__":
    main()
