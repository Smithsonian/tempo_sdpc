#! /usr/bin/env python3

import os, sys
import time
from datetime import date
import dateutil.parser
import sqlite3
import re
import argparse

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3, "defer":100}
DryRun = False
DB_Path = None

class Tokenizer:
    def __init__ (self):
        self.regex = re.compile (r'\s+=\s+')

    def tokens (self, line):
        s = line.rstrip(';\n')
        return self.regex.split(s)

def connect_database (mode):
    conn = sqlite3.connect ("file:{}?mode={}".format(DB_Path, mode), uri=True)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def get_product_table_names (cur):
    cur.execute ("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'");
    table_names = [item for t in cur.fetchall() for item in t]
    # Return table names in descending level order.
    # As a side-effect, this imposes an upload sequence, with L0 and RAW last
    table_names.sort(key=lambda x:x.split('_')[-1], reverse=True)
    return table_names

def table_files_matching_status (cur, table_name, asdc_status, limit=0, order='asc'):

    qual = "order by istart {}".format(order)

    if limit > 0:
        qual += " limit {}".format(limit)

    nc_query  = "select path from {table_name} where asdc_status == {asdc_status} {qual}".format(**locals())
    met_query = "select path from {table_name} where asdc_status_met == {asdc_status} {qual}".format(**locals())

    cur.execute (nc_query)
    nc_paths = [item for t in cur.fetchall() for item in t]
    cur.execute (met_query)
    met_paths = [item + ".met" for t in cur.fetchall() for item in t]

    return sorted (nc_paths + met_paths)

def files_matching_status (cur, asdc_status, **kwargs):
    """
    Returns a dict type, e.g:
        list["RAD_L1"] = list of RAD_L1 files with the specified asdc_status
    """
    table_names = get_product_table_names (cur)

    paths = {}
    for tbl in table_names:
        if tbl == "RAD_L1a":
            continue
        paths[tbl] = table_files_matching_status (cur, tbl, asdc_status, **kwargs)

    return paths

def table_name_for_file (filename):
    tok = filename.split('_')
    return '{}_{}'.format(tok[1], tok[2])

def table_name_for_file_raw (filename):
    return "RAW"

def set_file_stats (cur, table_name, file_basename, st):
    size  = st.st_size
    mtime = int(st.st_mtime)
    sql = "update {table_name} set mtime={mtime},size={size} where filename=\"{file_basename}\"".format (**locals())
    if DryRun:
        print(sql)
    else:
        cur.execute(sql)

def update_file_status (cur, filename, asdc_status, status_time, update_stat=False, disposition=None):
    file_basename = os.path.basename (filename)
    ext_split = os.path.splitext(file_basename)
    if '.nc' == ext_split[1]:
        status_var_name = 'asdc_status'
        table_name = table_name_for_file (file_basename)
    elif '.met' == ext_split[1]:
        status_var_name = 'asdc_status_met'
        file_basename = os.path.basename(ext_split[0])
        table_name = table_name_for_file (file_basename)
    elif '.tar' == ext_split[1]:
        status_var_name = 'asdc_status'
        table_name = table_name_for_file_raw (file_basename)
    else:
        print ('*** update_file_status: unsupported extension: {}'.format(file_basename))
        return

    if update_stat:
        set_file_stats (cur, table_name, file_basename, os.stat(filename))

    if asdc_status == Asdc_Status["uploaded"]:
        status_time_var = "asdc_upload_time"
    elif asdc_status == Asdc_Status["problem"] or asdc_status == Asdc_Status["accepted"]:
        status_time_var = "asdc_ingest_time"
    else:
        status_time_var = None

    if disposition is not None:
        disposition_str = ",asdc_disposition=\"{}\"".format(disposition)
    else:
        disposition_str = ""

    if status_time_var is None:
        sql = "update {table_name} set {status_var_name}={asdc_status}{disposition_str} where filename=\"{file_basename}\"".format(**locals())
    else:
        sql = "update {table_name} set {status_var_name}={asdc_status},{status_time_var}={status_time}{disposition_str} where filename=\"{file_basename}\"".format(**locals())

    if DryRun:
        print(sql)
    else:
        cur.execute (sql)

def count_files_matching_status (asdc_status):
    with connect_database("ro") as conn:
        table_lists = files_matching_status (conn.cursor(), asdc_status)
    num_files=0
    for table in table_lists.keys():
        num_files += len(table_lists[table])
    print(num_files)

def print_files_matching_status (asdc_status, **kwargs):
    with connect_database("ro") as conn:
        table_lists = files_matching_status (conn.cursor(), asdc_status, **kwargs)
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
    parsed_t = dateutil.parser.isoparse(tok[1])
    entry["time_stamp"] = int(parsed_t.timestamp())
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
                asdc_status = Asdc_Status["accepted"]
            else:
                asdc_status = Asdc_Status["problem"]
                num_bad += 1
            update_file_status (cur, entry["basename"], asdc_status, entry["time_stamp"], disposition=entry["disposition"])

    return num_bad

def process_longpan_files (longpan_file_list):
    with connect_database("rw") as conn:
        cur = conn.cursor()
        for longpan_file in longpan_file_list:
            try:
                num_bad = process_longpan (cur, longpan_file)
                if num_bad > 0:
                    print ('{} has {} bad files'.format(longpan_file, num_bad))
            except BaseException as e:
                print('An exception occurred: {}'.format(e))
                print ("Error processing file: {}".format(longpan_file))

def set_file_status (status, file_list, update_stat):
    with open(file_list, "r") as fp:
        files = fp.readlines()
    files = [f.strip() for f in files]
    status_time = int(time.time())
    with connect_database("rw") as conn:
        cur = conn.cursor()
        for f in files:
            update_file_status (cur, f, Asdc_Status[status], status_time, update_stat=update_stat)

def print_query (cur, sql):
    cur.execute (sql)
    while True:
        row = cur.fetchone()
        if row is None:
            break
        print(','.join(map(str,row)))

def print_report (asdc_status_name, ymd):
    asdc_status = Asdc_Status[asdc_status_name]
    if ymd:
        times="datetime(asdc_upload_time,\"unixepoch\"),datetime(asdc_ingest_time,\"unixepoch\")"
    else:
        times="asdc_upload_time,asdc_ingest_time"
    other_columns="asdc_disposition,path"
    print ("#{times},{other_columns}".format(**locals()))
    with connect_database("ro") as conn:
        cur = conn.cursor()
        table_names = get_product_table_names (cur)
        for tbl in table_names:
            if tbl == "RAD_L1a":
                continue
            sql = "select {times},{other_columns} from {tbl} where asdc_status = {asdc_status} or asdc_status_met = {asdc_status} order by asdc_upload_time".format (**locals())
            print_query (cur, sql)

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC file upload status')
    parser.add_argument('--num', metavar='STATUS', default=None,
                        help="Count files matching status: {}".format(Asdc_Status))
    parser.add_argument('--list', metavar='STATUS', default=None,
                        help="List files matching status: {}".format(Asdc_Status))
    parser.add_argument('--limit', metavar='LIMIT', default=0, type=int,
                        help="Maximum number of query results to list from each database table")
    parser.add_argument('--order', default='asc',
                        help="File path sort order (asc | desc)")
    parser.add_argument('--set', metavar=('STATUS','FILE_LIST',), default=None, nargs=2,
                        help="Set status of specified files")
    parser.add_argument('--pans', metavar='LONGPAN', default=None, nargs="*",
                        help="Process LONGPAN files (changes status 'uploaded' to 'accepted'|'problem')")
    parser.add_argument('--report', metavar='STATUS', default=None,
                        help="Print disposition report for files matching status")
    parser.add_argument('--ymd', action='store_true',
                        help="Print report times as YYYY-MM-DD HH:MM:SS")
    parser.add_argument('--stat', action='store_true',
                        help="Update file size, mtime")
    parser.add_argument('--dryrun', action='store_true',
                        help="Show actions, but don't modify the database")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    global DryRun
    DryRun = args.dryrun

    global DB_Path
    DB_Path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if DB_Path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    if args.num:
        count_files_matching_status (Asdc_Status[args.num])
    elif args.list:
        print_files_matching_status (Asdc_Status[args.list], limit=args.limit, order=args.order)
    elif args.set:
        set_file_status (args.set[0], args.set[1], args.stat)
    elif args.pans:
        process_longpan_files(args.pans)
    elif args.report:
        print_report (args.report, args.ymd)

if __name__ == "__main__":
    main()
