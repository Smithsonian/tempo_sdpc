#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import sqlite3
from datetime import datetime
import dateutil.parser
import dateutil.relativedelta
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import argparse

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def default_db_path():
    ancillary_root = os.getenv ("SDPC_ANCILLARY_ROOT")
    if ancillary_root == None:
        eprint ('*** Error: SDPC_ANCILLARY_ROOT is not set')
        sys.exit(1)
    return os.path.join (ancillary_root, "var/iers/iers.sqlite")

def get_file_keys (filename):
    nc = NetCDFFile (filename, "r")
    keys = {}
    keys["time_coverage_start"] = nc.getncattr ('time_coverage_start')
    nc.close()
    return keys

def time_info (t):
    tt = t.date().timetuple
    year = tt().tm_year
    yday = tt().tm_yday
    daytag = 1000 * year + yday
    info = {'year': year, 'yday': yday, 'daytag': daytag}
    return info

def run_sql_select (c, field, subst):
    cmd = "select path from 'IERS' where ({field} > :beg and {field} < :end) order by abs({field} - :tx)".format(**locals())
    c.execute(cmd, subst)
    rows = c.fetchall()
    rows = [r for r in rows if None not in r]
    return rows

def select_iers_file (c, window_days, keys):
    utc_string= keys["time_coverage_start"]
    t = dateutil.parser.isoparse (utc_string)
    dt = dateutil.relativedelta.relativedelta (days=window_days)
    dt1 = dateutil.relativedelta.relativedelta (days=1.0)

    # Don't select the IERS bulletin that exactly matches the target date!
    prev = time_info (t-dt1)
    beg  = time_info (t-dt)

    subst = {'tx':prev["daytag"], 'beg':beg["daytag"], 'end':prev["daytag"]}
    rows = run_sql_select (c, 'daytag', subst)

    return rows

def main():
    parser = argparse.ArgumentParser(description='Select an appropriate IERS bulletin A file')
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument ('--window', metavar='DAYS', default=60, type=int,
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

    if args.dbfile is None:
        db_path = default_db_path()
    else:
        db_path = args.dbfile

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    with sqlite3.connect (db_path) as conn:
        conn.execute("pragma foreign_keys=on")
        #conn.set_trace_callback(print)
        c = conn.cursor()
        iers = select_iers_file (c, args.window, file_keys)

    if len(iers) > 0:
        best_match = iers[0]
        print(best_match[0])

if __name__ == "__main__":
    main()
