#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import re
import time
import textwrap
import sqlite3
from datetime import date
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import yaml
import argparse

Sqlite_Trace = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_rad_info (path):
    with NetCDFFile (path, "r") as nc:
        tstart = nc.getncattr ('time_coverage_start_since_epoch')
        tend = nc.getncattr ('time_coverage_end_since_epoch')
        scan_num = nc.getncattr ('scan_num')

    info = {}
    info["tmid"] = 0.5*(tstart+tend)
    info["scan_num"] = scan_num
    return info

Rad_Basename_Parser = re.compile ("TEMPO_RAD_L1\w?_(?:|NRT_)V\d{2}_(\d{8}T\d{6}Z)_S\d{3}G\d{2}.nc")

def lookup_radiance_path (c, basename):

    g = Rad_Basename_Parser.match (basename)
    if g is None:
        return None

    timestamp = g.group(1)
    sql = "select path from 'RAD_L1' where filename like '%{}%';".format(timestamp)
    c.execute (sql)
    result = c.fetchone()
    if result is None:
        return None
    else:
        return result[0]

def select_matching_radref (c, window_days, window_hours, minwidth, thisscan, info):
    sec_per_day = 86400.0
    sec_per_hour = 3600.0
    tx = info["tmid"]
    htx = tx % sec_per_day
    tmax = window_days * sec_per_day + window_hours * sec_per_hour;

    if thisscan:
        scan_num_label = "%%S%03d%%" % (info["scan_num"])
        subst = {'beg':tx - 3*sec_per_hour, 'end':tx, 'scan_num_label': scan_num_label}
        cmd = 'select path from "RADREF_L1" where 0.5*(tstart + tend) between :beg and :end and filename like :scan_num_label order by tstart desc limit 1;'
    else:
        t1 = tx-tmax
        t2 = tx-sec_per_day
        subst = {'beg':min(t1,t2), 'end':max(t1,t2), 'minwidth':minwidth,
                 'tx':tx, 'htx':htx, 'window_days':window_days, 'window_hours':window_hours}
        cmd = """
              with cte as (
              select path,
                     abs(:tx - tstart)/86400.0 as day_diff,
                     abs(:htx - tstart%86400)/3600.0 as hr_diff
              from RADREF_L1
              where tstart between :beg and :end and num_mirror_pos >= :minwidth
              )
              select path from cte
              where day_diff < :window_days and hr_diff < :window_hours
              order by day_diff, hr_diff asc limit 1
              """
        cmd = textwrap.dedent(cmd)

    row = c.execute(cmd, subst)
    if row is None:
        path = ""
    else:
        path = row.fetchone()[0]

    return path

def connect_database (db_path):
    """
    For back-compatibility sqlite has foreign keys turned off by default,
    and foreign_keys=off is ALWAYS stored in the database, regardless of
    the runtime setting when the database was created.  For this reason,
    we apparently need to turn it on explicitly, each time the database
    connection is established.
    """
    conn = sqlite3.connect ("file:{}?mode=ro".format(db_path), uri=True, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    global Sqlite_Trace
    if Sqlite_Trace:
        conn.set_trace_callback(print)
    return conn

def config_file_defaults (filename):
    install_root = os.getenv ("SDPC_ROOT")
    if install_root == None:
        eprint ("*** Error: SDPC_ROOT is not set")
        sys.exit(1)

    pipe_root = os.getenv ("SDPC_PIPE_DIR")
    if pipe_root == None:
        eprint ("*** Error: SDPC_PIPE_DIR is not set")
        sys.exit(1)

    search_dirs = [pipe_root, install_root]

    cfg = None
    for dir in search_dirs:
        path = os.path.join (dir, filename)
        if os.path.isfile(path):
           try:
              with open(path) as fp:
                  cfg = yaml.load (fp, Loader=yaml.BaseLoader)
           except Exception as e:
               print('ERROR: {}'.format(e))

    return cfg

def main():

    cfg = config_file_defaults ("etc/select_radref.yml")
    if cfg is not None:
        days  = cfg["max_days_previous"]
        hours = cfg["max_time_offset"]
        minwidth = cfg["min_scan_width"]
    else:
        days = 2.0
        hours= 2.0
        minwidth = 600

    parser = argparse.ArgumentParser(description='Select the appropriate radiance reference file')
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument ('--days',default=days, type=float,
                         help="Acceptable offset in days")
    parser.add_argument ('--hours', default=hours, type=float,
                         help="Acceptable offset in hours")
    parser.add_argument ('--minwidth', default=minwidth, type=int,
                         help="Minimum acceptable scan width")
    parser.add_argument ('--thisscan', action='store_true', help="scan specific search")
    parser.add_argument ('--trace', action='store_true', help="trace sqlite query")
    parser.add_argument ('radiance_file', help="Basename of Level 1 radiance file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    global Sqlite_Trace
    Sqlite_Trace = args.trace

    if args.dbfile is None:
        db_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
        if db_path == None:
            eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
            sys.exit(1)
    else:
        db_path = args.dbfile

    rad_basename = os.path.basename (args.radiance_file)
    if "_NRT_" in rad_basename:
        rad_basename = rad_basename.replace ("_NRT", "")

    with connect_database (db_path) as conn:
        radiance_path = lookup_radiance_path (conn.cursor(), rad_basename)

    if radiance_path is None:
        print ('*** radiance file basename not in database: {}'.format (rad_basename))
        sys.exit(1)

    rad_info = get_rad_info (radiance_path)

    with connect_database (db_path) as conn:
        c = conn.cursor()
        try:
            radref = select_matching_radref (c, args.days, args.hours, args.minwidth, args.thisscan, rad_info)
        except:
            radref = ""

    print(radref)

if __name__ == "__main__":
    main()
