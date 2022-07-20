#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import math
import signal
from threading import Event
from subprocess import check_output

import sqlite3
from datetime import datetime
import dateutil.parser
import dateutil.relativedelta
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import argparse

Epoch = None
Epoch_Timet = None
Silent = False
DryRun = False

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

def emit_single_telemetry_only (t0, t1, destdir):
    global Epoch
    global Epoch_Timet
    timestamp = check_output (["ttime", "--epoch", Epoch, "--timestamp", "{}".format(t0)])
    timestamp = timestamp.decode("ascii")
    filename = 'TEMPO_INR_L1_V01_%s.nc' % (timestamp)
    path = os.path.join (destdir, filename)
    t0_struct = dateutil.parser.isoparse(timestamp)
    t0_timet = datetime.timestamp(t0_struct)
    if not Silent:
        print ('{}: telemetry-only -> {}'.format (datetime.now().isoformat(), path), flush=True)
    if not DryRun:
        # Atomically create file at 'path':
        temp_path = os.path.join (destdir, '.' + filename)
        with open (temp_path, "w") as fp:
            print ("{},{},{},{}".format(t0, t1, Epoch_Timet, t0_timet), file=fp)
        os.rename (temp_path, path)

def emit_n_telemetry_only (t_prev_iru, t_i, file_i, destdir, sig, wait_time):
    file_duration = 300.0     # seconds
    min_file_duration = 10.0

    with NetCDFFile (file_i, 'r') as nc:
        t_end = nc.time_coverage_end_since_epoch

    basename = os.path.basename(file_i)
    is_radiance = (basename.find('TEMPO_RAD') >= 0)

    if not is_radiance:
        t_i = t_end

    if t_i - t_prev_iru < file_duration:
        if is_radiance:
            if t_i - t_prev_iru > min_file_duration:
                emit_single_telemetry_only (t_prev_iru, t_i, destdir)
                sig.wait (wait_time)
            return t_end
        else:
            return t_prev_iru

    total_duration = t_i - t_prev_iru
    num_files = math.ceil(total_duration/ file_duration)
    mean_duration = total_duration / num_files

    for i in range(num_files):
        t0 = t_prev_iru + i * mean_duration
        t1 = t0 + mean_duration
        emit_single_telemetry_only (t0, t1, destdir)
        sig.wait (wait_time)

    return t_end

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
    parser.add_argument('--telemonly', action='store_true',
                        help="Trigger generation of telemetry-only radiance files")
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

    global Silent
    Silent = args.silent

    global DryRun
    DryRun = args.dryrun

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

    global Epoch
    with NetCDFFile (files[0], 'r') as nc:
        Epoch = nc.time_reference
    epoch_struct = dateutil.parser.isoparse(Epoch)
    global Epoch_Timet
    Epoch_Timet = datetime.timestamp(epoch_struct)

    t_prev = t[0]
    t_prev_iru = t_prev
    for i in range(len(files)):
        dt = t[i] - t_prev

        if args.wait is None:
            wait_time = dt
        else:
            wait_time = args.wait

        if sig.caught():
            print('\nCaught signal:  Resume using --start {}'.format(int(t[i])))
            break

        if args.telemonly:
            t_prev_iru = emit_n_telemetry_only (t_prev_iru, t[i], files[i], args.destdir, sig, wait_time)

        t_prev = t[i]
        src = files[i]
        dst = os.path.join (args.destdir, os.path.basename(src))
        if not Silent:
            print('{}: {} -> {}'.format(datetime.now().isoformat(), src, dst), flush=True)
        if not DryRun:
            os.symlink (src, dst)
        sig.wait (wait_time)

if __name__ == "__main__":
    main()
