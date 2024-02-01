#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import glob
import re
import time
from datetime import datetime
import argparse
import numpy as np

Verbose = False
Dryrun = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# We could read a more precise timestamp by opening the file,
# but parsing filename is sufficient, and is probably faster,
# since we only need one glob per directory.
def filename_start_time (base):
    # parse filename to get timestamp to whole seconds precision
    filename_regex = 'OR_ABI-L2-CMIPF-M\d{1}C\d{2}_G\d{2}_s(\d{13})\d_'
    fields = re.search (filename_regex, base)
    if fields is None:
        eprint ("*** Error: regex mismatch: {}".format(base))
        return None
    tstamp = fields.group(1)
    tstamp_obj = datetime.strptime(tstamp, '%Y%j%H%M%S')
    timet = time.mktime(tstamp_obj.timetuple())
    return timet

def utc_time_string (t):
    return datetime.utcfromtimestamp(t).strftime('%Y-%m-%dT%H:%M:%SZ')

def find_gaps (files, dt):

    times = []
    for f in files:
        times.append (filename_start_time (os.path.basename(f)))

    sorted_times = np.asarray(sorted(times))
    deltas = np.diff(sorted_times, 1)

    gap_indices = np.argwhere (deltas > dt)
    if len(gap_indices) == 0:
        return

    # print a CSV format line for each gap detected
    for k in gap_indices:
        i = k[0]
        t1 = sorted_times[i]
        t2 = sorted_times[i+1]
        print("{},{},{}".format(utc_time_string(t1), utc_time_string(t2), deltas[i]))

def main():
    default_delta = 1800.0
    parser = argparse.ArgumentParser(description='Find gaps in GOES CMI archive time coverage')
    parser.add_argument('--delta', type=float, default=default_delta, 
                        help="Minimum coverage gap [sec, default={}]".format (default_delta))
    parser.add_argument('dirs', help="Directory list", nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    files = []
    for dir in args.dirs:
        files += glob.glob(os.path.join (dir, "OR_ABI-L2-CMIPF-*.nc"))

    find_gaps (files, args.delta)

if __name__ == "__main__":
    main()
