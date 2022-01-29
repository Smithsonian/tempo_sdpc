#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import sqlite3
import argparse

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

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

def define_common_fields (fields):
    fields["path"] = "text"
    fields["daytag"] = "integer not null"
    fields["year"] = "integer not null"
    fields["yday"] = "integer not null"

def init_ims_table (table_name):
    fields = {}
    define_common_fields (fields)
    quals = "primary key(daytag), unique(daytag)"
    return Table_Type(table_name, fields, quals)

def insert_ims_entry (conn, table_name, entry):
    c = conn.cursor()
    ims = init_ims_table(table_name)
    ims.create(c)
    try:
        ims.new_entry (c, entry.keys(), entry.values())
        conn.commit()
        return 0
    except sqlite3.IntegrityError:
        eprint ('ERROR: duplicate primary key: daytag={}'.format(entry["daytag"]))
        return -1

def ims_entry_exists (conn, table_name, entry):
    c = conn.cursor()
    ims = init_ims_table(table_name)
    ims.create(c)
    c.execute ("select path from {} where daytag == {}".format (table_name, entry["daytag"]));
    path = c.fetchone()
    return path != None

def update_ims_entry (conn, table_name, entry):
    c = conn.cursor()
    sql = "update {} set path=\"{}\" where daytag={}".format (table_name, entry["path"], entry["daytag"]);
    try:
        c.execute (sql)
        conn.commit()
        return 0
    except:
        eprint ('ERROR: updating primary key: daytag={}'.format(entry["daytag"]))
        return -1

def process_file (conn, path):
    basename = os.path.basename (path)

    # example basename: ims2022027_1km_v1.3.nc.gz
    tok = basename.split('_')
    daytag = int(tok[0].strip('ims'))

    keys = {}
    keys["path"] = os.path.abspath(path)
    keys["daytag"] = daytag
    keys["year"] = int (daytag / 1000)
    keys["yday"] = daytag % 1000

    table_name = 'IMS'

    if ims_entry_exists (conn, table_name, keys):
        status = update_ims_entry (conn, table_name, keys)
    else:
        status = insert_ims_entry (conn, table_name, keys)

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
    return conn

def register_files (db_path, filenames):
    with connect_database (db_path) as conn:
        for fn in filenames:
            status = process_file (conn, fn)
            if status != 0:
                eprint('Error processing file: {}'.format(fn))

def main():
    parser = argparse.ArgumentParser(description='Register IMS files')
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if args.dbfile == None:
        eprint ('*** Error: DBFILE is not set')
        sys.exit(1)

    register_files (args.dbfile, args.filenames)

if __name__ == "__main__":
    main()
