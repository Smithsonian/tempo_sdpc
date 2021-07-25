#! /usr/bin/env python3

import os, sys
import time
from datetime import date
import sqlite3
import re
import argparse

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "accepted":2}
DryRun = False

class Tokenizer:
    def __init__ (self):
        self.regex = re.compile (r'\s+=\s+')

    def tokens (self, line):
        s = line.rstrip(';\n')
        return self.regex.split(s)

def connect_database ():
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    conn = sqlite3.connect (db_file_path)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def get_product_table_names (cur):
    cur.execute ("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'");
    table_names = [item for t in cur.fetchall() for item in t]
    return table_names

def table_files_matching_status (cur, table_name, asdc_status):
    cur.execute ("select path from {} where asdc_status == {} order by path".format(table_name, asdc_status))
    nc_paths = [item for t in cur.fetchall() for item in t]
    cur.execute ("select path from {} where asdc_status_met == {} order by path".format(table_name, asdc_status))
    met_paths = [item + ".met" for t in cur.fetchall() for item in t]
    return sorted (nc_paths + met_paths)

def files_matching_status (cur, asdc_status):
    """
    Returns a dict type, e.g:
        list["RAD_L1"] = list of RAD_L1 files with the specified asdc_status
    """
    table_names = get_product_table_names (cur)

    paths = {}
    for tbl in table_names:
        if tbl == "RAD_L1a":
            continue
        paths[tbl] = table_files_matching_status (cur, tbl, asdc_status)

    return paths

def table_name_for_file (filename):
    tok = filename.split('_')
    return '{}_{}'.format(tok[1], tok[2])

def update_file_status (cur, filename, asdc_status):
    file_basename = os.path.basename (filename)
    ext_split = os.path.splitext(file_basename)
    if '.nc' == ext_split[1]:
        status_var_name = 'asdc_status'
        table_name = table_name_for_file (file_basename)
    elif '.met' == ext_split[1]:
        status_var_name = 'asdc_status_met'
        file_basename = os.path.basename(ext_split[0])
        table_name = table_name_for_file (file_basename)
    elif '.raw' == ext_split[1]:
        status_var_name = 'asdc_status'
        table_name = 'RAW'
    else:
        print ('*** update_file_status: unsupported extension: {}'.format(file_basename))
        return
    sql = "update {} set {}={} where filename=\"{}\"".format(table_name, status_var_name, asdc_status, file_basename)
    if DryRun:
        print(sql)
    else:
        cur.execute (sql)

def count_files_matching_status (asdc_status):
    with connect_database() as conn:
        table_lists = files_matching_status (conn.cursor(), asdc_status)
    num_files=0
    for table in table_lists.keys():
        num_files += len(table_lists[table])
    print(num_files)

def print_files_matching_status (asdc_status):
    with connect_database() as conn:
        table_lists = files_matching_status (conn.cursor(), asdc_status)
    for table in table_lists.keys():
        for f in table_lists[table]:
            print(f)

def longpan_header (thefile, parse):
    # MESSAGE_TYPE = LONGPAN
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'MESSAGE_TYPE' or tok[1] != 'LONGPAN':
        print ('*** invalid file header: {}'.format(line))
        return -1
    # NO_OF_FILES = $num_files
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'NO_OF_FILES':
        print ('*** invalid file header: {}'.format(line))
        return -1
    return int(tok[1])

def longpan_entry (thefile, parse):
    entry = {}
    # FILE_DIRECTORY = $asdc_dir_path
    line = thefile.readline()
    tok = parse.tokens (line)
    # FILE_NAME = $basename
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'FILE_NAME':
        print ('*** invalid entry: {}'.format(line))
        return None
    entry["basename"] = tok[1]
    # DISPOSITION = $disposition
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'DISPOSITION':
        print ('*** invalid entry: {}'.format(line))
        return None
    entry["disposition"] = tok[1].strip('" ')
    # TIME_STAMP = $time_stamp
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'TIME_STAMP':
        print ('*** invalid entry: {}'.format(line))
        return None
    entry["time_stamp"] = tok[1]
    return entry

def process_longpan(cur, longpan_file):
    """
    A LONGPAN file has a 2-line header:
        MESSAGE_TYPE = LONGPAN;
        NO_OF_FILES = $num_files;
    followed by $num_files sets of 4-line entries:
        FILE_DIRECTORY = $asdc_dir_path;
        FILE_NAME = $basename;
        DISPOSITION = $disposition;
        TIME_STAMP = $time_stamp;
    where $time_stamp has the form YYYY-MM-DDTHH:MM:SSZ

    On success, $disposition = "SUCCESSFUL", otherwise there was a problem
    """
    parse = Tokenizer()
    num_bad = 0
    with open (longpan_file, "r") as thefile:
        num_files = longpan_header (thefile, parse)
        if num_files < 0:
            print ('skipping file: {}'.format(longpan_file))
            return
        for i in range(num_files):
            entry = longpan_entry (thefile, parse)
            if entry == None:
                break
            if entry["disposition"] == "SUCCESSFUL":
                update_file_status (cur, entry["basename"], Asdc_Status["accepted"])
            else:
                update_file_status (cur, entry["basename"], Asdc_Status["problem"])
                num_bad += 1

    return num_bad

def process_longpan_files (longpan_file_list):
    with connect_database() as conn:
        cur = conn.cursor()
        for longpan_file in longpan_file_list:
            try:
                num_bad = process_longpan (cur, longpan_file)
                if num_bad > 0:
                    print ('{} has {} bad files'.format(longpan_file, num_bad))
            except:
                print ("Error processing file: {}".format(longpan_file))

def set_file_status (status, file_list):
    with open(file_list, "r") as fp:
        files = fp.readlines()
    files = [f.strip() for f in files]
    with connect_database() as conn:
        cur = conn.cursor()
        for f in files:
            update_file_status (cur, f, Asdc_Status[status])

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC file upload status')
    parser.add_argument('--num', metavar='STATUS', default=None,
                        help="Count files matching status: {}".format(Asdc_Status))
    parser.add_argument('--list', metavar='STATUS', default=None,
                        help="List files matching status: {}".format(Asdc_Status))
    parser.add_argument('--set', metavar=('STATUS','FILE_LIST',), default=None, nargs=2,
                        help="Set status of specified files")
    parser.add_argument('--pans', metavar='LONGPAN', default=None, nargs="*",
                        help="Process LONGPAN files (changes status 'pending' to 'accepted'|'problem')")
    parser.add_argument('--dryrun', action='store_true',
                        help="Show actions, but don't modify the database")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    global DryRun
    DryRun = args.dryrun

    if args.num:
        count_files_matching_status (Asdc_Status[args.num])
    elif args.list:
        print_files_matching_status (Asdc_Status[args.list])
    elif args.set:
        set_file_status (args.set[0], args.set[1])
    elif args.pans:
        process_longpan_files(args.pans)

if __name__ == "__main__":
    main()
