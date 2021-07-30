#! /usr/bin/env python3

import os, sys
import time
import argparse
import csv

def keep_path (filter, hour, path):
    basename = os.path.basename (path)
    pattern = filter['pattern']
    begin_exclude_hour = filter['begin_exclude_hour']
    end_exclude_hour = filter['end_exclude_hour']
    n = len(pattern)
    for i in range(n):
        if basename.find(pattern[i]) < 0:
            continue
        try:
            beg_hour = float(begin_exclude_hour[i])
        except:
            beg_hour = 0
        try:
            end_hour = float(end_exclude_hour[i])
        except:
            end_hour = 24.0
        if beg_hour <= hour and hour <= end_hour:
            return False
    return True

def strip_comments(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

# strip comments and leading/trailing spaces
def read_csv_columns (path):
    with open (path,'r') as fp:
        reader = csv.reader (strip_comments(fp))
        headers = next (reader, None)
        headers = [h.strip() for h in headers]
        column = {}
        for h in headers:
            column[h] = []
        for row in reader:
            for h, v in zip(headers, row):
                column[h].append(v.strip())
    return column

def main():
    parser = argparse.ArgumentParser(description='Filter a list of files')
    parser.add_argument('--filter', help="Filter definition file (CSV)")
    parser.add_argument('infile', metavar='FILELIST', help="Input file list")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if args.filter == None:
        parser.print_usage(sys.stderr)
        sys.exit(1)

    filter = read_csv_columns (args.filter)

    # min resolution is good enough
    t = time.localtime()
    hour = t.tm_hour + t.tm_min/60.0

    orig_file = args.infile + '.orig'
    os.rename (args.infile, orig_file)

    with open (orig_file, 'r') as fp:
        input_paths = fp.readlines()

    with open (args.infile, 'w') as fp:
        for path in input_paths:
            if keep_path (filter, hour, path):
                fp.write(path)

if __name__ == "__main__":
    main()
