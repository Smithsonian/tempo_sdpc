#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import re
import os, sys
import sqlite3
import argparse

DryRun = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Tokenizer:
    def __init__ (self):
        self.regex = re.compile (r'\s+=\s+')

    def tokens (self, line):
        s = line.rstrip(';\n')
        return self.regex.split(s)

class Table_Type:

    def __init__ (self, table_name, fields, quals):
        self.table_name = table_name
        self.field_defs = fields
        field_list = ','.join('{} {}'.format (k,fields[k]) for k in fields.keys())
        self.create_cmd = "CREATE TABLE IF NOT EXISTS {table_name} ({field_list}, {quals});".format (**locals())

    def create(self, cur):
        """Create an empty table"""
        cur.execute (self.create_cmd)

    def new_entry(self, cur, names, values):
        """Insert a table entry"""
        self.create(cur)
        name_list=','.join(names)
        value_slots = ','.join(['?']*len(values))
        cmd = "INSERT INTO {} ({}) VALUES ({})".format(self.table_name, name_list, value_slots)
        value_string_tuple = tuple(str(v) for v in values)
        cur.execute (cmd, value_string_tuple)

def init_file_table (table_name):
    fields = {}
    fields["rowid"] = "integer"
    fields["asdc_status"] = "integer"
    fields["filename"] = "text"
    fields["path"] = "text"
    quals = "primary key(rowid)"
    return Table_Type(table_name, fields, quals)

def insert_file_entry (conn, table_name, entry):
    c = conn.cursor()
    tbl = init_file_table(table_name)
    tbl.create(c)
    tbl.new_entry (c, entry.keys(), entry.values())
    conn.commit()
    return 0

def file_entry_exists (conn, table_name, entry):
    c = conn.cursor()
    tbl = init_file_table(table_name)
    tbl.create(c)
    c.execute ("select path from {} where filename == \"{}\"".format (table_name, entry["filename"]));
    path = c.fetchone()
    return path != None

def update_file_entry (conn, table_name, entry):
    c = conn.cursor()
    sql = "update {} set path=\"{}\" where filename=\"{}\"".format (table_name, entry["path"], entry["filename"]);
    c.execute (sql)
    conn.commit()
    return 0

def process_file (conn, table_name, path):
    basename = os.path.basename (path)

    keys = {}
    keys["asdc_status"] = 0
    keys["filename"] = basename
    keys["path"] = os.path.abspath(path)

    if file_entry_exists (conn, table_name, keys):
        status = update_file_entry (conn, table_name, keys)
    else:
        status = insert_file_entry (conn, table_name, keys)

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
    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def files_matching_status1 (cur, table_name, asdc_status):
    cur.execute ("select path from {} where asdc_status == {} order by path".format(table_name, asdc_status))
    paths = [item for t in cur.fetchall() for item in t]
    return sorted (paths)

def update_file_status (cur, filename, table_name, asdc_status):
    basename = os.path.basename(filename)
    sql = "update {} set asdc_status={} where filename=\"{}\"".format(table_name, asdc_status, basename)
    if DryRun:
        print(sql)
    else:
        cur.execute (sql)

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

def process_longpan(cur, table_name, longpan_file, status_dict):
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
                update_file_status (cur, entry["basename"], table_name, status_dict["accepted"])
            else:
                update_file_status (cur, entry["basename"], table_name, status_dict["problem"])
                num_bad += 1

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

    def files_matching_status (self, asdc_status):
        with connect_database (self.db_path) as conn:
            file_list = files_matching_status1 (conn.cursor(), self.table_name, asdc_status)
        return file_list

    def set_file_status (self, status, file_list):
        with open(file_list, "r") as fp:
            files = fp.readlines()
        files = [f.strip() for f in files]
        with connect_database(self.db_path) as conn:
            cur = conn.cursor()
            for f in files:
                update_file_status (cur, f, self.table_name, self.status_dict[status])

    def process_longpan_files (self, longpan_file_list):
        with connect_database(self.db_path) as conn:
            cur = conn.cursor()
            for longpan_file in longpan_file_list:
                try:
                    num_bad = process_longpan (cur, self.table_name, longpan_file, self.status_dict)
                    if num_bad > 0:
                        print ('{} has {} bad files'.format(longpan_file, num_bad))
                except:
                    print ("Error processing file: {}".format(longpan_file))

def main():
    status_dict = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "accepted":2}

    parser = argparse.ArgumentParser(description='Track ASDC upload status of ancillary data files using {}'.format(status_dict.keys()))
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument('--num', metavar='STATUS', default=None,
                        help="Count files matching status")
    parser.add_argument('--list', metavar='STATUS', default=None,
                        help="List files matching status:")
    parser.add_argument('--set', metavar=('STATUS','FILE_LIST',), default=None, nargs=2,
                        help="Set status of specified files")
    parser.add_argument('--dryrun', action='store_true',
                        help="Show actions, but don't modify the database")
    parser.add_argument('--add', metavar='FILE', default=None, nargs="*",
                        help="Add new files to the database")
    parser.add_argument('--pans', metavar='LONGPAN', default=None, nargs="*",
                        help="Process LONGPAN files (changes status 'pending' to 'accepted'|'problem')")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if args.dbfile == None:
        eprint ('*** Error: DBFILE is not set')
        sys.exit(1)

    global DryRun
    DryRun = args.dryrun

    db = Db_File_Type (args.dbfile, 'File_Table', status_dict)

    if args.num:
        file_list = db.files_matching_status (status_dict[args.num])
        print(len(file_list))
    elif args.list:
        file_list = db.files_matching_status (status_dict[args.list])
        for f in file_list:
            print(f)
    elif args.set:
        db.set_file_status (args.set[0], args.set[1])
    elif args.add:
        db.register_files (args.add)
    elif args.pans:
        db.process_longpan_files(args.pans)

if __name__ == "__main__":
    main()
