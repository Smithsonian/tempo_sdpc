#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
from datetime import datetime, timezone
import argparse
import numpy as np
import sqlite3

Verbose = False
Dryrun = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def utc_time_string (t):
    return datetime.fromtimestamp(t).strftime('%Y-%m-%dT%H:%M:%SZ')

def connect_database (path):
    return sqlite3.connect ("file:{}?mode=ro".format(path), uri=True, timeout=20.0)

def print_coverage_gaps (dbfile, delta):
    sql = \
"""
    with cte as (
  select tstart,
          lead(tstart, 1, tstart) over (order by tstart) as tnext,
          lead(tstart, 1, tstart) over (order by tstart) - tstart as diff
  from File_Table order by tstart
)
select tstart, tnext, diff
from cte
where diff > {}
""".format(delta)

    with connect_database (dbfile) as conn:
        cur = conn.execute (sql)
        rows = cur.fetchall()

    print("# file: {}:".format(dbfile))
    print("# num_gaps: {}".format(len(rows)))
    print("begin_utc,end_utc,delta_sec")
    for row in rows:
        (t1,t2,dt) = row
        print("{},{},{}".format(utc_time_string(t1), utc_time_string(t2), dt))

def main():
    default_delta = 1800.0
    parser = argparse.ArgumentParser(description='Find gaps in GOES CMI archive time coverage')
    parser.add_argument('--delta', type=float, default=default_delta,
                        help="Minimum coverage gap [sec, default={}]".format (default_delta))
    parser.add_argument('dbfiles', metavar='DBFILES',
                        help="sqlite dbfile list", nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    print("# Created: {}".format(datetime.now()))

    for dbfile in args.dbfiles:
        print_coverage_gaps (dbfile, args.delta)

if __name__ == "__main__":
    main()
