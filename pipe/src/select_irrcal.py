#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import sqlite3
from netCDF4 import Dataset as NetCDFFile
import argparse

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_db_path():
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE_L1")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE_L1 is not set')
        sys.exit(1)

    return db_file_path

def get_tstart (filename):
    nc = NetCDFFile (filename, "r")
    tstart = nc.getncattr ('time_coverage_start_since_epoch')
    nc.close()
    return tstart

def table_exists (cur, table_name):
    cur.execute ("SELECT name FROM sqlite_master WHERE type='table' AND name='{}';".format(table_name))
    result = cur.fetchone()
    return result != None

def select_irrcal (cur, table_name, tstart):
    istart = int(tstart)
    cur.execute("select path from '{table_name}' where istart = {istart}".format (**locals()))
    result = cur.fetchone()
    if result is None:
        path = ""
    else:
        path = result[0]
    return path

def main():
    parser = argparse.ArgumentParser(description='Select the matching Level 2 irradiance calibration file')
    parser.add_argument ('-m', '--molecule', metavar="MOLECULE", default=None, help="Molecule abbreviation (case sensitive)")
    parser.add_argument ('irr_file', help="Path to Level 1 irradiance file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if args.molecule is None:
        parser.print_usage(sys.stderr)
        sys.exit(1)

    tstart = get_tstart (args.irr_file)

    db_path = get_db_path()

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    with sqlite3.connect ("file:{}?mode=ro".format(db_path), uri=True, timeout=20.0) as conn:
        conn.execute("pragma foreign_keys=on")
        #conn.set_trace_callback(print)
        cur = conn.cursor()
        table_name = "IRR{}_L2".format(args.molecule)
        if table_exists (cur, table_name):
            path = select_irrcal (cur, table_name, tstart)
        else:
            path = ""

    print(path)

if __name__ == "__main__":
    main()
