#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import math
import csv
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

def connect_database (dbfile):
    conn = sqlite3.connect ("file:{}?mode=ro".format(dbfile), uri=True, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def time_filter (beg, end):
    if beg is None and end is None:
        where = ""
    elif end is None:
        where = "where {beg} <= time_coverage_start_since_epoch".format (**locals())
    elif beg is None:
        where = "where time_coverage_end_since_epoch < {end}".format (**locals())
    else:
        where= "where {beg} <= time_coverage_start_since_epoch and time_coverage_end_since_epoch < {end}".format (**locals())
    return where

def tables_in_db (dbfile, want_tables):
    with connect_database(dbfile) as conn:
        cur = conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        existing_tables = cur.fetchall()
    result = []
    for t in existing_tables:
        if t[0] in want_tables:
            result.append (t[0])
    if len(result) == 0:
        result = None
    return result

def tables_union_query (tables, beg, end):
    where = time_filter (beg, end)
    sql_parts = []
    for tbl in tables:
        sql = "select time_coverage_start_since_epoch,path from '{tbl}' {where}".format (**locals())
        sql_parts.append (sql)
    if len(sql_parts) > 1:
        sql = " union ".join(sql_parts)
    return sql + ' order by time_coverage_start_since_epoch'

def select_files (dbfile, want_level1, beg, end):

    if want_level1:
        tables = ["RAD_L1"]
    else:
        image_types = ['DRK', 'DRKL', 'IRR', 'IRRL', 'IRRR', 'RAD', 'RADT']
        tables = [s + '_L0' for s in image_types]

    tables = tables_in_db (dbfile, tables)
    if tables is None:
        return None

    sql = tables_union_query (tables, beg, end)
    if sql is None:
        return None

    with connect_database (dbfile) as conn:
        cur = conn.execute(sql)
        rows = cur.fetchall()

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

  def handler(self, signum, frame):
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

def read_csv_file_list (csvpath):
    if not os.path.isfile(csvpath):
        return None
    tstart = []
    path = []
    with open (csvpath, 'r') as fp:
        reader = csv.reader(fp)
        for row in reader:
            tstart.append (int(row[0]))
            path.append (row[1])

    tstart = np.asarray(tstart)
    path = np.asarray(path)
    indices = np.argsort(tstart)
    return {"tstart":tstart[indices], "path":path[indices]}

class TrackDB:

    dbfile = None

    def __create_if_not_exists (self):
        table_name = "File_Table"
        fields = {}
        fields["istart"] = "integer not null"
        fields["filename"] = "text"
        fields["path"] = "text"
        fields["tstamp"] = "integer"
        quals = "primary key(istart)"
        field_list = ','.join ('{} {}'.format (k, fields[k]) for k in fields.keys())
        sql = "create table if not exists {table_name} ({field_list}, {quals});".format (**locals())
        # journal_mode=WAL mode is persistent
        with sqlite3.connect (self.dbfile) as conn:
            conn.execute("pragma journal_mode=WAL")
            conn.execute (sql)

    def __define_db_path (self, dbfile):
        if '/' in dbfile:
            self.dbfile = dbfile
            return
        pipe_dbfile = os.getenv ("SDPC_ARCHIVE_DBFILE")
        if pipe_dbfile is None:
            raise Exception("SDPC_ARCHIVE_DBFILE is not set")
        dir = os.path.dirname (pipe_dbfile)
        if not os.path.isdir(dir):
            raise Exception("Nonexistent directory: {}".format(dir))
        self.dbfile = os.path.join (dir, dbfile)
        if self.dbfile == pipe_dbfile:
            raise Exception("Tracking database file must differ from SDPC_ARCHIVE_DBFILE")

    def __init__ (self, dbfile):
        self.__define_db_path (dbfile)
        self.__create_if_not_exists ()

    def has_file (self, path):
        filename = os.path.basename(path)
        sql = "select exists(select 1 from File_Table where filename = '{filename}')".format(**locals())
        with sqlite3.connect (self.dbfile) as conn:
            cur = conn.cursor()
            cur.execute (sql)
            result = cur.fetchone()
        return int(result[0]) != 0

    def insert_file (self, tstart, path):
        filename = os.path.basename(path)
        istart = int(tstart)
        tstamp = int(time.time())
        sql = "insert into File_Table (istart,filename,path,tstamp) values ({istart},'{filename}','{path}',{tstamp})".format(**locals())
        with sqlite3.connect (self.dbfile) as conn:
            conn.execute(sql)

def main():
    parser = argparse.ArgumentParser(description='Deliver time-ordered Level 0 or Level 1 image data to a target directory')
    parser.add_argument('--dbfile', metavar='DBFILE',
                        help="Path to sqlite database containing source files")
    parser.add_argument('--csvfile', metavar='CSVFILE',
                        help="Path to CSV file containing source file list (overrides dbfile)")
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
    parser.add_argument('--notrack', action='store_true',
                        help="Disable tracking of source files")
    parser.add_argument('--trackdbfile', default="reprocess_input_files.sqlite",
                        help="Specify non-default filename for sqlite source file tracking database")
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

    file_dict = None
    if args.csvfile is not None:
        file_dict = read_csv_file_list (args.csvfile)
    elif args.dbfile is not None:
        file_dict = select_files (args.dbfile, args.level1, args.start, args.end)

    if file_dict is None:
        eprint ('*** Error: No source files to deliver: dbfile={} csvfile={}'.format(args.dbfile, args.csvfile))
        sys.exit(1)

    files = file_dict['path']
    times = np.asarray(file_dict['tstart'])

    # Check the files actually exist
    num_files = 0
    for f in files:
        if not os.path.isfile(f):
            eprint ('*** Cannot access file: {}'.format(f))
            sys.exit(1)
        num_files += 1

    print('Files input: {} (all exist)'.format(num_files), flush=True)

    if args.notrack:
        trackdb = None
        do_tracking = False
    else:
        trackdb = TrackDB (args.trackdbfile)
        do_tracking = True

    sig = Signal_Catcher()

    global Epoch
    with NetCDFFile (files[0], 'r') as nc:
        Epoch = nc.time_reference
    epoch_struct = dateutil.parser.isoparse(Epoch)
    global Epoch_Timet
    Epoch_Timet = datetime.timestamp(epoch_struct)

    num_linked = 0
    num_skipped = 0

    t_prev = times[0]
    t_prev_iru = t_prev
    for i in range(len(files)):
        dt = times[i] - t_prev

        if args.wait is None:
            wait_time = dt
        else:
            wait_time = args.wait

        if sig.caught():
            print('\nCaught signal')
            if args.dbfile is not None:
                print('Resume using --start {}'.format(int(times[i])))
            break

        if args.telemonly:
            t_prev_iru = emit_n_telemetry_only (t_prev_iru, times[i], files[i], args.destdir, sig, wait_time)

        t_prev = times[i]
        src = files[i]
        dst = os.path.join (args.destdir, os.path.basename(src))

        if do_tracking:
            if trackdb.has_file (src):
                eprint ('already processed: {}'.format(src), flush=True)
                num_skipped += 1
                continue

        if not Silent:
            print('{}: {} -> {}'.format(datetime.now().isoformat(), src, dst), flush=True)

        if not DryRun:
            if do_tracking:
                trackdb.insert_file (times[i], src)
            os.symlink (src, dst)
            num_linked += 1

        sig.wait (wait_time)

    if not Silent:
         print ('Files processed: total:{} linked:{} skipped:{}'.format(num_files, num_linked, num_skipped), flush=True)

if __name__ == "__main__":
    main()
