#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import glob
import re
import time
from datetime import datetime, timezone
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
    # help strptime by adding a timezone specification
    tstamp = fields.group(1)+'+0000'
    tstamp_obj = datetime.strptime(tstamp, '%Y%j%H%M%S%z')
    timet = time.mktime(tstamp_obj.timetuple())
    return timet

def utc_time_string (t):
    return datetime.fromtimestamp(t).strftime('%Y-%m-%dT%H:%M:%SZ')

def find_gaps (files, dt):

    times = []
    for f in files:
        times.append (filename_start_time (os.path.basename(f)))

    sorted_times = sorted(np.asarray(times))
    deltas = np.diff(sorted_times, 1)

    gap_indices = np.argwhere (deltas > dt)
    if len(gap_indices) == 0:
        return;

    t1 = []
    t2 = []
    dt = []
    for k in gap_indices:
        i = k[0]
        t1.append(sorted_times[i])
        t2.append(sorted_times[i+1])
        dt.append(deltas[i])

    gaps = {"beg":t1, "end":t2, "delta": dt};

    return gaps

def main():
    default_delta = 1800.0
    parser = argparse.ArgumentParser(description='Find gaps in GOES CMI archive daily time coverage')
    parser.add_argument('--delta', type=float, default=default_delta,
                        help="Minimum coverage gap [sec, default={}]".format (default_delta))
    parser.add_argument('dirs', help="Directory list", nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    for dir in args.dirs:
        files = glob.glob(os.path.join (dir, "OR_ABI-L2-CMIPF-*.nc"))
        gaps = find_gaps (files, args.delta)
        if gaps is not None:
            t1 = [utc_time_string(s) for s in gaps["beg"]]
            t2 = [utc_time_string(s) for s in gaps["end"]]
            dt = gaps["delta"]
            n = len(dt)
            print("gap: {}".format(dir))
            for i in range(n):
                print("{},{},{}".format(t1[i], t2[i], dt[i]))

if __name__ == "__main__":
    main()
