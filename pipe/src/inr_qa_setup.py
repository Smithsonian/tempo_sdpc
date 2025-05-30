#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import re
import argparse
import sqlite3
from dateutil import parser as dateparser
import time

from netCDF4 import Dataset

TraceSQL = False

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, flush=True, **kwargs)

def connect_database (dbfile):
    conn = sqlite3.connect ("file:{}?mode=ro".format(dbfile), uri=True, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    if TraceSQL:
        conn.set_trace_callback(print)
    return conn

def get_radiance_paths (dbfile, t1, t2, scan_id=None):
    # Unix time_t at the TEMPO epoch (1980-01-06T00:00:00Z)
    tempo_epoch_timet = 315964800
    # TEMPO time window bounds
    istart1 = t1 - tempo_epoch_timet
    istart2 = t2 - tempo_epoch_timet

    # Optional filter on scan number
    if scan_id is None:
        scan_filter = ''
    else:
        scan_filter = "and scan_num = %d" % scan_id

    with connect_database (dbfile) as conn:
        cur = conn.cursor()
        cur.execute ("select path from RAD_L1 where istart between {} and {} {} order by istart".format(istart1, istart2, scan_filter))
        rows = cur.fetchall()

    paths = [t[0] for t in rows]
    print ("{} radiance files: {}".format(len(paths), dbfile))
    return paths

def get_goes_paths (dbfile, t1, t2):
    dbfile = os.path.expandvars (dbfile)
    with connect_database (dbfile) as conn:
        cur = conn.cursor()
        cur.execute ("select path from File_Table where tstart between {} and {}".format(t1, t2))
        rows = cur.fetchall()

    paths = [t[0] for t in rows]
    if len(paths) == 0:
        eprint ("Error: No GOES imagery in time range: {}-{}".format(t1, t2))
        return None
    else:
        print ("{} GOES files: {}".format (len(paths), dbfile))

    return paths

def classify_goes_file (path):
    with Dataset (path, "r") as nc:
        platform_id = nc.getncattr ('platform_ID')
    return int(platform_id.strip('G'))

def create_config_file_text (work_dir, radiance_dir, goes_dir):
    """
    The seemingly pointless first line MUST be present,
    and the following lines MUST be in the order specified:
        - Diary_path --> diary_directory
        - granule_path --> granule_directory
        - ref_im_path --> ref_im_directory
        - OutputParentFolder --> output_parent_directory
        - TempMatFiles --> temp_MAT_files_directory
    """
    config_file_lines = [\
    "variable_name=variable_path",
    "Diary_path={work_dir}/diary/",
    "granule_path={radiance_dir}/",
    "ref_im_path={goes_dir}/",
    "OutputParentFolder={work_dir}/output/",
    "TempMatFiles={work_dir}/tmp/"]
    config_file_text = "\n".join (config_file_lines) + "\n"
    return config_file_text.format(**locals())

def make_symlinks (dir, paths):
    for p in paths:
        dir_p = os.path.join (dir, os.path.basename(p))
        os.symlink (p, dir_p)

def initialize_working_directory (work_dir, rad_paths, goes_paths, goes_west_id, goes_east_id):
    # populate subdirectory with TEMPO radiances
    radiance_dir = os.path.join (work_dir, "radiances")
    os.mkdir (radiance_dir)
    make_symlinks (radiance_dir, rad_paths)
    # populate subdirectory with GOES imagery
    goes_dir = os.path.join (work_dir, "goes")
    os.mkdir (goes_dir)
    make_symlinks (goes_dir, goes_paths)
    # create working directories
    for d in ["diary", "output", "tmp"]:
        dir = os.path.join (work_dir, d)
        os.mkdir (dir)
    # generate config file
    config_file_path = os.path.join (work_dir, "config.txt")
    config_file_text = create_config_file_text (work_dir, radiance_dir, goes_dir)
    with open (config_file_path, "w") as fp:
        fp.write(config_file_text)
    with open (os.path.join (work_dir, "GOES_East"), "w") as fp:
        fp.write(format(goes_east_id))
    with open (os.path.join (work_dir, "GOES_West"), "w") as fp:
        fp.write(format(goes_west_id))
    print("initialized working directory: {}".format(work_dir))

def main():
    parser = argparse.ArgumentParser(description='Initialize working directory for INR quality check')
    parser.add_argument('--dbfile', help="TEMPO archive sqlite database path")
    parser.add_argument('--dir', help="Target directory path (must exist)")
    parser.add_argument('--scan', default=None, type=int, help="(optional) radiance scan number")
    parser.add_argument('date', help="Date of interest (YYYY-MM-DD)")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    dbfile = args.dbfile
    work_dir = args.dir
    scan_id = args.scan
    date = args.date

    if not os.path.isdir (work_dir):
        eprint ("Error: Nonexistent working directory: {}".format(work_dir))
        sys.exit(1)

    """
    Select input files for the specified date.
    The time window spans the maximum interval within instrument safety constraints.
    """
    # Unix time_t on the specified date, at approximately the earliest possible scan start:
    tm_tuple = dateparser.parse(date + 'T04:00:00-06:00').timetuple()

    # time_t window bounds [seconds]
    t1 = time.mktime(tm_tuple)
    t2 = t1 + 16*3600

    # Collect TEMPO radiances
    rad_paths = get_radiance_paths (dbfile, t1, t2, scan_id=scan_id)
    if rad_paths is None:
        sys.exit(1)

    # Collect GOES imagery
    goes_east = "$SDPC_ANCILLARY_ROOT/var/goes/cmieast.sqlite"
    goes_west = "$SDPC_ANCILLARY_ROOT/var/goes/cmiwest.sqlite"
    goes_east_paths = get_goes_paths (goes_east, t1, t2)
    if goes_east_paths is None:
        sys.exit(1)
    goes_east_id = classify_goes_file (goes_east_paths[0])

    goes_west_paths = get_goes_paths (goes_west, t1, t2)
    if goes_west_paths is None:
        sys.exit(1)
    goes_west_id = classify_goes_file (goes_west_paths[0])

    goes_paths = goes_east_paths + goes_west_paths

    initialize_working_directory (work_dir, rad_paths, goes_paths, goes_west_id, goes_east_id)

if __name__ == "__main__":
    main()
