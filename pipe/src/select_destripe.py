#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import time
import sqlite3
from datetime import date
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import yaml
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

def get_rad_info (path):
    with NetCDFFile (path, "r") as nc:
        tstart = nc.getncattr ('time_coverage_start_since_epoch')
        tend = nc.getncattr ('time_coverage_end_since_epoch')

    info = {}
    info["tmid"] = 0.5*(tstart+tend)
    return info

def lookup_radiance_path (c, basename):
    sql = "select path from 'RAD_L1' where filename='{basename}';".format(**locals())
    c.execute (sql)
    result = c.fetchone()
    if result is None:
        return None
    else:
        return result[0]

def select_matching_destripe (c, molecule, window_days, window_hours, minwidth, info):
    sec_per_day = 86400.0
    sec_per_hour = 3600.0
    tx = info["tmid"]
    htx = tx % sec_per_day
    tmax = window_days * sec_per_day + window_hours * sec_per_hour;

    # First, get the candidate list, with the newest listed first
    subst = {'beg':tx-tmax, 'end':tx-sec_per_day, 'minwidth':minwidth}
    cmd = "select path,0.5*(tstart+tend) from 'DSTR{}_L2' where :beg <= 0.5*(tstart+tend) and 0.5*(tstart+tend) <= :end and num_mirror_pos >= :minwidth order by tstart desc;".format(molecule)
    c.execute(cmd, subst)

    # Select from the list
    path = ""
    min_t_delta = window_days
    min_h_delta = window_hours / 24.0
    tol = min_t_delta * min_h_delta
    for row in c:
        tmid = row[1]
        t_delta = abs(tx - tmid) / sec_per_day
        h_delta = abs(htx - (tmid % sec_per_day)) / sec_per_day
        score = t_delta * h_delta
        if score < tol:
            path = row[0]
            tol = score

    return path

def connect_database (db_path):
    """
    For back-compatibility sqlite has foreign keys turned off by default,
    and foreign_keys=off is ALWAYS stored in the database, regardless of
    the runtime setting when the database was created.  For this reason,
    we apparently need to turn it on explicitly, each time the database
    connection is established.
    """
    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
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

    cfg = config_file_defaults ("etc/select_destripe.yml")
    if cfg is not None:
        days  = cfg["max_days_previous"]
        hours = cfg["max_time_offset"]
        minwidth = cfg["min_scan_width"]
    else:
        days = 2.0
        hours= 2.0
        minwidth = 600

    parser = argparse.ArgumentParser(description='Select the appropriate destriping correction file')
    parser.add_argument ('--molecule',default="HCHO",
                         help="Molecule symbol")
    parser.add_argument ('--days',default=days, type=float,
                         help="Acceptable offset in days")
    parser.add_argument ('--hours', default=hours, type=float,
                         help="Acceptable offset in hours")
    parser.add_argument ('--minwidth', default=minwidth, type=int,
                         help="Minimum acceptable scan width")
    parser.add_argument ('radiance_file', help="Basename of Level 1 radiance file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    db_path = get_db_path()

    with connect_database (db_path) as conn:
        radiance_path = lookup_radiance_path (conn.cursor(), args.radiance_file)

    if radiance_path is None:
        print ('*** radiance file not in database: {}'.format (args.radiance_file))
        sys.exit(1)

    rad_info = get_rad_info (radiance_path)

    with connect_database (db_path) as conn:
        c = conn.cursor()
        try:
            destripe = select_matching_destripe (c, args.molecule, args.days, args.hours, args.minwidth, rad_info)
        except:
            destripe = ""

    print(destripe)

if __name__ == "__main__":
    main()
