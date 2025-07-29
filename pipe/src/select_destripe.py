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

def get_radiance_start_time (c, basename):
    c.execute ("select istart from 'RAD_L1' where filename='{}';".format(basename))
    result = c.fetchone()
    if result is None:
        return None
    else:
        return result[0]

def find_most_recent_irradiance (c, istart):
    c.execute ("select istart from 'IRR_L1' where istart < {} order by istart desc limit 1".format(istart))
    result = c.fetchone()
    if result is None:
        return None
    else:
        return result[0]

def select_suitable_destripe (c, molecule, window_days, istart, istart_solar):
    """
    Look for a pre-computed destriping file generated since the most recent
    solar calibration observation.  If none exists, then an alternate destriping
    method will be used.
    """
    dt_max = window_days * 86400
    dt_irr = istart - istart_solar

    if dt_irr < dt_max:
        istart_min = istart_solar
    else:
        istart_min = istart - dt_max

    cmd = "select path from 'DSTR{}_L2' where istart > {} order by istart desc limit 1".format(molecule, istart_min)
    c.execute(cmd)
    path = c.fetchone()
    if path is None:
        return None
    else:
        return path[0]

def connect_database (db_path):
    conn = sqlite3.connect ("file:{}?mode=ro".format(db_path), uri=True, timeout=20.0)
    """
    For back-compatibility sqlite has foreign keys turned off by default,
    and foreign_keys=off is ALWAYS stored in the database, regardless of
    the runtime setting when the database was created.  For this reason,
    we apparently need to turn it on explicitly, each time the database
    connection is established.
    """
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
    else:
        days = 1.0

    parser = argparse.ArgumentParser(description='Select the appropriate destriping correction file')
    parser.add_argument ('--molecule', help="Molecule symbol")
    parser.add_argument ('--days',default=days, type=float, help="Acceptable offset in days")
    parser.add_argument ('radiance_path', help="Path to Level 1 radiance file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    db_path = get_db_path()

    with NetCDFFile (args.radiance_path, "r") as nc:
        istart = int(nc.time_coverage_start_since_epoch)

    with connect_database (db_path) as conn:
        istart_solar = find_most_recent_irradiance (conn.cursor(), istart)

    if istart_solar is None:
        print ('*** no irradiance file in database! (this should never happen)')
        sys.exit(1)

    with connect_database (db_path) as conn:
        try:
            destripe = select_suitable_destripe (conn.cursor(), args.molecule, args.days, istart, istart_solar)
        except:
            destripe = ""

    print(destripe)

if __name__ == "__main__":
    main()
