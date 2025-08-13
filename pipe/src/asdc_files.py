#! /usr/bin/env python3

# This program is intended to manage a database of files of a single type
# that are uploaded to ASDC.  It's primary purpose is to help track the
# ASDC upload status.  The files tracked may have custom fields that
# may be used in queries performed by other programs. This program
# can define the custom fields for such entries, but it doesn't
# use them in any way.

# for eprint definition
from __future__ import print_function

import re
import os, sys
import time
import calendar
import sqlite3
import argparse

DryRun = False
TraceSQL = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def abi_start_time_utc_tuple (base):
    # parse filename to get UTC time-tuple to whole seconds precision
    filename_regex = r'OR_ABI-L2-CMIPF-M\d{1}C\d{2}_G\d{2}_s(\d{13})\d_'
    fields = re.search (filename_regex, base)
    if fields is None:
        eprint ("*** Error: regex mismatch: {}".format(base))
        return None
    return time.strptime (fields.group(1), '%Y%j%H%M%S')

def abi_start_time_timet (base):
    """
    Compute unix time_t value for ABI filename's UTC start time
    """
    tp_utc = abi_start_time_utc_tuple (base)
    return calendar.timegm(tp_utc)

class Tokenizer:
    def __init__ (self):
        self.regex = re.compile (r'\s+=\s+')

    def tokens (self, line):
        s = line.rstrip(';\n')
        return self.regex.split(s)

class File_Type (object):
    def __init__ (self, regex, fields_method, entry_method):
        self.regex = re.compile(regex)
        self.fields = fields_method
        self.entry = entry_method

# Each file type that needs custom database fields
# should have an entry in this dict:
Filetype_Dict = {}

# Each Filetype_Dict entry should provide a regular
# expression to classify filenames, plus two functions:
#    *) one to define the custom fields,
#    *) one to populate the custom fields given the path
#       to a specific file

# IMS -----------------------------------
def ims_fields ():
    fields = {}
    fields["daytag"] = "integer not null"
    fields["year"] = "integer not null"
    fields["yday"] = "integer not null"
    return fields

def ims_entry (path):
    basename = os.path.basename(path)
    # example basename: ims2022027_1km_v1.3.nc.gz
    tok = basename.split('_')
    daytag = int(tok[0].strip('ims'))
    fields = {}
    fields["daytag"] = daytag
    fields["year"] = int (daytag / 1000)
    fields["yday"] = daytag % 1000
    return fields

def abi_fields ():
    fields = {}
    fields["tstart"] = "integer not null"
    return fields

def abi_entry (path):
    basename = os.path.basename(path)
    fields = {}
    fields["tstart"] = abi_start_time_timet (basename)
    return fields

Filetype_Dict["ims"] = File_Type(r"ims\d{7,7}_1km_v\d.\d.nc", ims_fields, ims_entry)
Filetype_Dict["abi"] = File_Type(r'OR_ABI-L2-CMIPF-M\d{1}C\d{2}_G\d{2}_s\d{13}\d_', abi_fields, abi_entry)

def classify_filename (path, *args, **kwargs):
    basename = os.path.basename(path)
    for key, value in Filetype_Dict.items():
        if value.regex.match(basename) is not None:
            return key
    return None

def file_fields (path, *args, **kwargs):
    type_string = classify_filename (path)
    if type_string is None:
        return {}
    else:
        return Filetype_Dict[type_string].fields(*args, **kwargs)

def file_entry (path, *args, **kwargs):
    type_string = classify_filename (path)
    if type_string is None:
        return {}
    else:
        return Filetype_Dict[type_string].entry(path, *args, **kwargs)

class Table_Type:
    def __init__ (self, table_name, fields, quals):
        self.table_name = table_name
        self.field_defs = fields
        field_list = ','.join('{} {}'.format (k,fields[k]) for k in fields.keys())
        self.create_cmd = "CREATE TABLE IF NOT EXISTS {table_name} ({field_list}, {quals});".format (**locals())
        self.create_trigger_cmd = \
        'create trigger if not exists update_status_time update of asdc_status on {table_name} '\
        'begin' \
        '  update {table_name} set asdc_status_time=current_timestamp where filename = old.filename; ' \
        'end;'.format (**locals())
        self.create_index_cmd = \
        "create unique index if not exists filename_index_{table_name} on {table_name}(filename);".format (**locals())

    def create(self, cur):
        """Create an empty table"""
        cur.execute (self.create_cmd)
        cur.execute (self.create_trigger_cmd)
        cur.execute (self.create_index_cmd)

    def new_entry(self, cur, names, values):
        """Insert a table entry"""
        self.create(cur)
        name_list=','.join(names)
        value_slots = ','.join(['?']*len(values))
        cmd = "INSERT INTO {} ({}) VALUES ({})".format(self.table_name, name_list, value_slots)
        value_string_tuple = tuple(str(v) for v in values)
        cur.execute (cmd, value_string_tuple)

def init_file_table (table_name, fields_for_file_type):
    fields = {}
    fields["rowid"] = "integer"
    fields["timestamp"] = "datetime default current_timestamp"
    fields["asdc_status_time"] = "datetime"
    fields["asdc_status"] = "integer"
    fields["filename"] = "text"
    fields["path"] = "text"
    quals = "primary key(rowid)"
    fields.update(fields_for_file_type)
    return Table_Type(table_name, fields, quals)

def insert_file_entry (cur, tbl, entry):
    tbl.new_entry (cur, entry.keys(), entry.values())
    return 0

def file_entry_exists (cur, table_name, entry):
    cur.execute ("select path from {} where filename == \"{}\"".format (table_name, entry["filename"]));
    path = cur.fetchone()
    return path != None

def update_file_entry (cur, table_name, entry):
    sql = "update {} set path=\"{}\" where filename=\"{}\"".format (table_name, entry["path"], entry["filename"]);
    cur.execute (sql)
    return 0

def process_file (conn, table_name, path):
    basename = os.path.basename (path)

    fields = {}
    fields["asdc_status"] = 0
    fields["filename"] = basename
    fields["path"] = os.path.abspath(path)
    fields.update (file_entry (path))

    cur = conn.cursor()
    tbl = init_file_table (table_name, file_fields (path))
    tbl.create (cur)

    if file_entry_exists (cur, table_name, fields):
        status = update_file_entry (cur, table_name, fields)
    else:
        status = insert_file_entry (cur, tbl, fields)

    conn.commit()

    if status < 0:
        eprint('ERROR: processing file {}'.format(path))

    return status

def connect_database (db_path):
    """
    For back-compatibility sqlite has foreign keys turned off by default,
    and foreign_keys=off is ALWAYS stored in the database, regardless of
    the runtime setting when the database was created.  For this reason,
    we apparently need to turn it on explicitly, each time the database
    connection is established.
    """
    conn = sqlite3.connect (db_path, timeout=20.0)
    conn.execute("pragma foreign_keys=on")
    if TraceSQL:
        conn.set_trace_callback(print)
    return conn

def files_matching_status1 (cur, table_name, asdc_status, limit=0):
    if limit > 0:
        query = "select path from {} where asdc_status == {} order by path limit {}".format(table_name, asdc_status, limit)
    else:
        query = "select path from {} where asdc_status == {} order by path".format(table_name, asdc_status)
    cur.execute (query)
    paths = [item for t in cur.fetchall() for item in t]
    return sorted (paths)

def update_file_status (conn, table_name, params):
    sql = "update {} set asdc_status=:asdc_status where filename=:filename".format(table_name)
    if DryRun:
        print(sql)
    else:
        conn.executemany (sql, params)

def query_file_status (cur, filename, table_name, missing):
    basename = os.path.basename(filename)
    query = "select asdc_status from {} where filename=\"{}\"".format (table_name, basename)
    cur.execute (query)
    asdc_status = cur.fetchone()
    if asdc_status == None:
        return missing
    else:
        return asdc_status[0]

def longpan_header (thefile, parse):
    # MESSAGE_TYPE = LONGPAN
    line = thefile.readline()
    tok = parse.tokens (line)
    if tok[0] != 'MESSAGE_TYPE' or tok[1].strip('" ') != 'LONGPAN':
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
    entry["basename"] = tok[1].strip('" ')
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

def process_longpan(db_path, table_name, longpan_file, status_dict):
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

    params = []
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
                asdc_status = status_dict["accepted"]
            else:
                asdc_status = status_dict["problem"]
                num_bad += 1
            params.append({"asdc_status":asdc_status, "filename":entry["basename"]})

    with connect_database(db_path) as conn:
        update_file_status (conn, table_name, params)

    return num_bad

class Db_File_Type:

    def __init__ (self, db_path, table_name, status_dict):
        self.db_path = db_path
        self.table_name = table_name
        self.status_dict = status_dict

    def register_files (self, filenames):
        with connect_database (self.db_path) as conn:
            for fn in filenames:
                status = process_file (conn, self.table_name, fn)
                if status != 0:
                    eprint('Error processing file: {}'.format(fn))

    def files_matching_status (self, asdc_status, **kwargs):
        with connect_database (self.db_path) as conn:
            file_list = files_matching_status1 (conn.cursor(), self.table_name, asdc_status, **kwargs)
        return file_list

    def set_file_status (self, status, file_list):
        with open(file_list, "r") as fp:
            files = fp.readlines()
        files = [f.strip() for f in files]
        asdc_status = self.status_dict[status]
        params = []
        for f in files:
            params.append({"asdc_status":asdc_status, "filename":os.path.basename(f)})
        with connect_database(self.db_path) as conn:
            update_file_status (conn, self.table_name, params)

    def print_file_status (self, filenames):
        with connect_database (self.db_path) as conn:
            cur = conn.cursor()
            for fn in filenames:
                asdc_status = query_file_status (cur, fn, self.table_name, self.status_dict["unknown"])
                print(asdc_status)

    def process_longpan_files (self, longpan_file_list):
        for longpan_file in longpan_file_list:
            try:
                num_bad = process_longpan (self.db_path, self.table_name, longpan_file, self.status_dict)
                if num_bad > 0:
                    print ('{} has {} bad files'.format(longpan_file, num_bad))
            except:
                print ("Error processing file: {}".format(longpan_file))

def main():
    status_dict = {"unknown": -3, "nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3}

    parser = argparse.ArgumentParser(description='Track ASDC upload status of ancillary data files using {}'.format(status_dict.keys()))
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument('--num', metavar='STATUS', default=None,
                        help="Count files matching status")
    parser.add_argument('--list', metavar='STATUS', default=None,
                        help="List files matching status:")
    parser.add_argument('--limit', metavar='LIMIT', default=0, type=int,
                        help="Maximum number of query results to list")
    parser.add_argument('--set', metavar=('STATUS','FILE_LIST',), default=None, nargs=2,
                        help="Set status of specified files")
    parser.add_argument('--dryrun', action='store_true',
                        help="Show actions, but don't modify the database")
    parser.add_argument('--trace', action='store_true',
                        help="Trace SQL transactions")
    parser.add_argument('--add', metavar='FILE', default=None, nargs="*",
                        help="Add new files to the database")
    parser.add_argument('--pans', metavar='LONGPAN', default=None, nargs="*",
                        help="Process LONGPAN files (changes status 'pending' to 'accepted'|'problem')")
    parser.add_argument('--status', metavar='FILE', default=None, nargs="*",
                        help="Query status of specified files")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if args.dbfile == None:
        eprint ('*** Error: DBFILE is not set')
        sys.exit(1)

    global TraceSQL
    TraceSQL = args.trace

    global DryRun
    DryRun = args.dryrun

    db = Db_File_Type (args.dbfile, 'File_Table', status_dict)

    if args.num:
        file_list = db.files_matching_status (status_dict[args.num])
        print(len(file_list))
    elif args.list:
        file_list = db.files_matching_status (status_dict[args.list], limit=args.limit)
        for f in file_list:
            print(f)
    elif args.set:
        db.set_file_status (args.set[0], args.set[1])
    elif args.add:
        db.register_files (args.add)
    elif args.pans:
        db.process_longpan_files(args.pans)
    elif args.status:
        db.print_file_status (args.status)

if __name__ == "__main__":
    main()
