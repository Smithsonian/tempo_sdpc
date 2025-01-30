#! /usr/bin/env python3

import os, sys
import time
from datetime import date
import dateutil.parser
import sqlite3
import re
import fnmatch
import argparse
import traceback

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3, "defer":100}
DryRun = False
DB_Path = None
TraceSQL = False

PDR_DB_Path = None  # used only to process SHORTPAN files

Uploads_Excluded = ["RAD_L1a", "RADT_L1a"]

class Tokenizer:
    def __init__ (self):
        self.regex = re.compile (r'\s+=\s+')

    def tokens (self, line):
        s = line.rstrip(';\n')
        return self.regex.split(s)

def connect_database (mode):
    conn = sqlite3.connect ("file:{}?mode={}".format(DB_Path, mode), uri=True, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    if TraceSQL:
        conn.set_trace_callback(print)
    return conn

def __connect_database (mode, dbfile):
    conn = sqlite3.connect ("file:{}?mode={}".format(dbfile, mode), uri=True, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    if TraceSQL:
        conn.set_trace_callback(print)
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
        if tbl in Uploads_Excluded:
            continue
        paths[tbl] = table_files_matching_status (cur, tbl, asdc_status, **kwargs)

    return paths

def table_name_for_file (path_or_basename):
    basename = os.path.basename (path_or_basename)
    ext_split = os.path.splitext (basename)
    if '.tar' == ext_split[1]:
        return "RAW"
    tok = basename.split('_')
    return '{}_{}'.format(tok[1], tok[2])

def product_file_path (nc_or_met_file_path):
    ext_split = os.path.splitext (nc_or_met_file_path)
    if '.met' == ext_split[1]:
        return ext_split[0]
    elif nc_or_met_file_path.endswith ('.cmr.json'):
        return nc_or_met_file_path.replace ('.cmr.json', '')
    else:
        return nc_or_met_file_path

def update_file_status (cur, table_name, filename, asdc_status, status_time, update_stat=False, disposition=None):
    fields = {}

    file_basename = os.path.basename (filename)
    ext_split = os.path.splitext(file_basename)
    if '.nc' == ext_split[1]:
        fields["asdc_status"] = asdc_status
    elif '.met' == ext_split[1]:
        fields["asdc_status_met"] = asdc_status
        file_basename = ext_split[0]
    elif file_basename.endswith ('.cmr.json'):
        fields["asdc_status_met"] = asdc_status
        file_basename = file_basename.replace ('.cmr.json','')
    elif '.tar' == ext_split[1]:
        fields["asdc_status"] = asdc_status
    else:
        print ('*** update_file_status: unsupported extension: {}'.format(file_basename))
        return

    if update_stat:
        path = product_file_path (filename)
        if os.path.isfile (path):
            st = os.stat(path)
            fields["mtime"] = int(st.st_mtime)
            fields["size"] = st.st_size

    if asdc_status == Asdc_Status["uploaded"]:
        fields["asdc_upload_time"] = status_time
    elif asdc_status == Asdc_Status["problem"] or asdc_status == Asdc_Status["accepted"]:
        fields["asdc_ingest_time"] = status_time

    if disposition is not None:
        fields["asdc_disposition"] = disposition

    set_named_fields = ",".join([k + "=:" + k for k in fields.keys()])

    sql = "update {table_name} set {set_named_fields} where filename=\"{file_basename}\"".format (**locals())

    if DryRun:
        print(sql)
    else:
        cur.execute (sql, fields)

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

def shortpan_entry (thefile, parse):
    entry = {}
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'DISPOSITION':
        print ('*** invalid file header: {}'.format(line))
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

def longpan_header (thefile, parse):
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

def get_product_filenames_from_pdr_file (pdr_path):
    with open (pdr_path, "r") as fp:
        lines = fp.readlines()
    lines = [s.strip() for s in lines]
    lines = list(filter (None, lines))
    # Extract the product filenames from the PDR file
    p = re.compile (r"\s*FILE_ID\s*=\s*(..*);$")
    product_files = []
    for s in lines:
        m = re.match (p, s)
        if m is not None:
            product_files.append(m.group(1))
    return product_files

def process_shortpan (cur, thefile, parse, pan_file):
    entry = shortpan_entry (thefile, parse)
    if entry is None:
        print ('skipping file: {}'.format(pan_file))
        return
    # Can we look up the corresponding PDR file to get product filenames?
    global PDR_DB_Path
    if PDR_DB_Path is None:
        print ('skipping SHORTPAN: {} (pdrdbfile=None)'.format(pan_file))
        return
    # Replace PAN file extension (.pan or .PAN) to get the PDR file's local path
    pdr_file = os.path.basename(re.sub (".pan", ".PDR", pan_file, flags=re.IGNORECASE))
    pdr_query = "select path from File_Table where filename == '{}'".format(pdr_file)
    with __connect_database ("ro", PDR_DB_Path) as conn:
        pdr_cur = conn.cursor()
        pdr_cur.execute (pdr_query)
        result = pdr_cur.fetchone()
        if result is None:
            print ('SQL query yields empty result: {}'.format(pdr_query))
            return
        else:
            pdr_path = result[0]
    # Try to read the local PDR file
    if not os.path.isfile (pdr_path):
        print ('skipping SHORTPAN: {} (cannot open corresponding PDR file: {})'.format(pan_file, pdr_path))
        return
    product_files = get_product_filenames_from_pdr_file (pdr_path)
    if entry["disposition"] == "SUCCESSFUL":
        asdc_status = Asdc_Status["accepted"]
    else:
        asdc_status = Asdc_Status["problem"]
        print ('unsuccessful shortpan {} references {} files'.format(pan_file, len(product_files)))
    # Assign the shortpan disposition to the remaining product files from the PDR.
    for basename in product_files:
        table_name = table_name_for_file (basename)
        update_file_status (cur, table_name, basename, asdc_status, entry["time_stamp"], disposition=entry["disposition"])

def process_longpan (cur, thefile, parse, pan_file):
    num_files = longpan_header (thefile, parse)
    if num_files < 0:
        print ('skipping file: {}'.format(pan_file))
        return

    num_bad = 0
    for i in range(num_files):
        entry = longpan_entry (thefile, parse)
        if entry == None:
            break
        if entry["disposition"] == "SUCCESSFUL":
            asdc_status = Asdc_Status["accepted"]
        else:
            asdc_status = Asdc_Status["problem"]
            num_bad += 1
        table_name = table_name_for_file (entry["basename"])
        update_file_status (cur, table_name, entry["basename"], asdc_status, entry["time_stamp"], disposition=entry["disposition"])

    if num_bad > 0:
        print ('{} has {} bad files'.format(pan_file, num_bad))
    return

def process_pan (cur, pan_file):
    """
    The generic SIPS ICD, 423-41-57_ICD_ECSSIPS_RevK_Final.pdf, says:
       There are two forms of the PAN, a short (Table 3-14) and a long
       (Table 3-15) form. The short form of the PAN is sent to
       acknowledge that all files have been successfully transferred,
       or to report errors that are not specific to individual files
       but which have precluded processing of any and all files (e.g.,
       transfer failure). If all files in a request do not have the
       same disposition, the long form of this message is employed.

    A SHORTPAN file just has 3 lines:
        MESSAGE_TYPE = SHORTPAN;
        DISPOSITION = $disposition;
        TIME_STAMP = $time_stamp

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
    with open (pan_file, "r") as thefile:
        line = thefile.readline()
        tok = parse.tokens (line)
        if tok[0] != 'MESSAGE_TYPE':
            print ('*** unexpected file header: {}'.format(line))
            return -1
        message_type = tok[1].strip('"')  # <-- Because people can't be bothered follow specifications
        if  message_type == 'LONGPAN':
            process_longpan (cur, thefile, parse, pan_file)
        elif message_type == 'SHORTPAN':
            process_shortpan (cur, thefile, parse, pan_file)
        else:
            print ('*** unexpected file header: {}'.format(line))
            return -1

    return 0

def process_pan_files (pan_file_list):
    for pan_file in pan_file_list:
        try:
            with connect_database("rw") as conn:
                process_pan (conn.cursor(), pan_file)
        except BaseException as e:
            traceback.print_tb(e.__traceback__)
            print('An exception occurred: {}'.format(e))
            print("Error processing file: {}".format(pan_file))

def set_file_status (status, file_list, update_stat):
    with open(file_list, "r") as fp:
        files = fp.readlines()
    files = [f.strip() for f in files]
    # filter out empty strings
    files = list(filter (None, files))
    status_time = int(time.time())

    # To streamline transactions and minimize the duration each connection is held open,
    # we process the files by table, starting a new database connection for each table.

    table_lists = {}
    for f in files:
        table_name = table_name_for_file (f)
        if table_name in table_lists.keys():
            table_lists[table_name].append(f)
        else:
            table_lists[table_name] = [f]

    for table_name in table_lists.keys():
        with connect_database("rw") as conn:
            for f in table_lists[table_name]:
                update_file_status (conn.cursor(), table_name, f, Asdc_Status[status], status_time, update_stat=update_stat)

def print_query (cur, sql):
    cur.execute (sql)
    for row in cur:
        print(','.join(map(str,row)))

def print_report (asdc_status_name, ymd, limit=0):
    asdc_status = Asdc_Status[asdc_status_name]
    if ymd:
        times="datetime(asdc_upload_time,\"unixepoch\"),datetime(asdc_ingest_time,\"unixepoch\")"
    else:
        times="asdc_upload_time,asdc_ingest_time"
    other_columns="asdc_disposition,path"
    limit_qual = ""
    if limit > 0:
        limit_qual = "limit {}".format(limit)
    print ("#{times},{other_columns}".format(**locals()))
    with connect_database("ro") as conn:
        cur = conn.cursor()
        table_names = get_product_table_names (cur)
        for tbl in table_names:
            if tbl in Uploads_Excluded:
                continue
            sql = "select {times},{other_columns} from {tbl} where asdc_status = {asdc_status} or asdc_status_met = {asdc_status} order by asdc_upload_time {limit_qual}".format (**locals())
            print_query (cur, sql)

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC file upload status')
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
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
    parser.add_argument('--pans', metavar='PAN', default=None, nargs="*",
                        help="Process PAN files (changes status 'uploaded' to 'accepted'|'problem')")
    parser.add_argument('--pdrdbfile', metavar='PDRDBFILE', default=None,
                        help="sqlite database path used to track PAN files")
    parser.add_argument('--report', metavar='STATUS', default=None,
                        help="Print disposition report for files matching status")
    parser.add_argument('--ymd', action='store_true',
                        help="Print report times as YYYY-MM-DD HH:MM:SS")
    parser.add_argument('--stat', action='store_true',
                        help="Update file size, mtime")
    parser.add_argument('--dryrun', action='store_true',
                        help="Show actions, but don't modify the database")
    parser.add_argument('--trace', action='store_true',
                        help="Trace SQL statements")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    global TraceSQL
    TraceSQL = args.trace

    global DryRun
    DryRun = args.dryrun

    global DB_Path
    if args.dbfile is None:
        DB_Path = os.getenv ("SDPC_ARCHIVE_DBFILE")
        if DB_Path == None:
            eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
            sys.exit(1)
    else:
        DB_Path = args.dbfile

    global PDR_DB_Path
    PDR_DB_Path = args.pdrdbfile

    if args.num:
        count_files_matching_status (Asdc_Status[args.num])
    elif args.list:
        print_files_matching_status (Asdc_Status[args.list], limit=args.limit, order=args.order)
    elif args.set:
        set_file_status (args.set[0], args.set[1], args.stat)
    elif args.pans:
        process_pan_files(args.pans)
    elif args.report:
        print_report (args.report, args.ymd, limit=args.limit)

if __name__ == "__main__":
    main()
