#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import signal
from threading import Event

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

def time_filter (beg, end):
    if beg is None and end is None:
        where = ""
    elif end is None:
        where = "where time_coverage_start_since_epoch >= {beg}".format (**locals())
    elif beg is None:
        where = "where time_coverage_end_since_epoch < {end}".format (**locals())
    else:
        where= "where time_coverage_start_since_epoch >= {beg} and time_coverage_end_since_epoch < {end}".format (**locals())
    return where

def level1_query (c, beg, end):
    where = time_filter (beg, end)
    sql = "select time_coverage_start_since_epoch,path from 'RAD_L1' {where} order by time_coverage_start_since_epoch".format (**locals())
    return sql

def level0_tables(c):
    image_types = ['DRK', 'DRKL', 'IRR', 'IRRL', 'IRRR', 'RAD', 'RADT']
    tables = [s + '_L0' for s in image_types]
    c.execute("SELECT name FROM sqlite_master WHERE type='table'")
    rows = c.fetchall()
    result = []
    for r in rows:
        if r[0] in tables:
            result.append (r[0])
    if len(result) == 0:
        result = None
    return result

def level0_query (c, beg, end):
    tables = level0_tables (c)
    if tables is None:
        return None
    where = time_filter (beg, end)
    sql_parts = []
    for tbl in tables:
        sql = "select time_coverage_start_since_epoch,path from '{tbl}' {where}".format (**locals())
        sql_parts.append (sql)
    sql = " union ".join(sql_parts)
    return sql + ' order by time_coverage_start_since_epoch'

def select_files (c, want_level1, beg, end):

    if want_level1:
        sql = level1_query (c, beg, end)
    else:
        sql = level0_query (c, beg, end)

    if sql is None:
        return None

    c.execute(sql)
    rows = c.fetchall()
    tstart = []
    path = []
    for r in rows:
        tstart.append(r[0])
        path.append(r[1])

    return {'tstart':tstart, 'path':path}

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGTERM, self.handler)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

def main():
    parser = argparse.ArgumentParser(description='Deliver time-ordered Level 0 or Level 1 image data to a target directory')
    parser.add_argument('--dbfile', metavar='DBFILE',
                        help="Sqlite database path")
    parser.add_argument('--wait', default=None, type=float,
                        help="Time interval [sec] between files")
    parser.add_argument('--level1', action='store_true',
                        help="Feed Level 1 radiance data")
    parser.add_argument('--start', metavar='STARTTIME', default=None,
                        help="Start time [TAI sec since TEMPO epoch]")
    parser.add_argument('--end', metavar='ENDTIME', default=None,
                        help="End time [TAI sec since TEMPO epoch]")
    parser.add_argument('--silent', action='store_true',
                        help="Minimize output")
    parser.add_argument('--dryrun', action='store_true',
                        help="Dry run for testing")
    parser.add_argument('destdir', metavar='DESTDIR',
                        help="Destination directory to receive symbolic links")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    dbfile = args.dbfile

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    with sqlite3.connect (dbfile) as conn:
        conn.execute("pragma foreign_keys=on")
        #conn.set_trace_callback(print)
        c = conn.cursor()
        file_dict = select_files (c, args.level1, args.start, args.end)

    if file_dict is None:
        eprint ('No files selected')
        return

    files = file_dict['path']
    t = np.asarray(file_dict['tstart'])

    sig = Signal_Catcher()

    t_prev = t[0]
    for i in range(len(files)):
        if args.wait is None:
            dt = t[i] - t_prev
        else:
            dt = args.wait
        sig.wait(dt)
        if sig.caught():
            print('\nCaught signal:  Resume using --start {}'.format(int(t[i])))
            break
        t_prev = t[i]
        src = files[i]
        dst = os.path.join (args.destdir, os.path.basename(src))
        if not args.silent:
            print('{}: {} -> {}'.format(datetime.now().isoformat(), src, dst), flush=True)
        if not args.dryrun:
            os.symlink (src, dst)

if __name__ == "__main__":
    main()
