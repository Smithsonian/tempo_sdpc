#! /usr/bin/python

import os, sys
import time
import sqlite3
from datetime import date
from netCDF4 import Dataset as NetCDFFile

Radiance_Files = ["RAD_L1a", "RAD_L1b"]
Radiance_Products = ["CLDRR", "HCHO", "NO2", "O3T", "O3P"]
Radiance_Derived_Files = [s + "_L2" for s in Radiance_Products] \
                       + [s + "_L3" for s in Radiance_Products]

Radiance_File_Attributes = ["time_coverage_start_since_epoch", "time_coverage_end_since_epoch",
                            "scan_num", "scan_type", "granule_num"]

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

def init_radiance_table (table_name):
    fields = {}
    fields["istart"] = "integer not null"
    fields["scan_type"] = "integer not null"
    fields["scan_num"] = "integer not null"
    fields["granule_num"] = "integer not null"
    fields["time_coverage_start_since_epoch"] = "float not null"
    fields["time_coverage_end_since_epoch"] = "float not null"
    fields["filename"] = "text"
    fields["mtime"] = "float"
    fields["size"] = "integer"
    fields["mirror_pos_beg"] = "integer not null"
    fields["mirror_pos_end"] = "integer not null"
    quals = "primary key(istart), unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_radiance_product_table (table_name):
    fields = {}
    fields["istart"] = "integer not null"
    fields["mtime"] = "float"
    fields["size"] = "integer"
    fields["filename"] = "text"
    quals = "unique(istart), foreign key (istart) references {}(istart)".format('RAD_L1b')
    return Table_Type(table_name, fields, quals)

def init_other_product_table (table_name):
    fields = {}
    fields["istart"] = "integer not null"
    fields["mtime"] = "float"
    fields["size"] = "integer"
    fields["filename"] = "text"
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def insert_radiance_entry (conn, table_name, entry):
    with conn:
        c = conn.cursor()
        rad = init_radiance_table(table_name)
        rad.create(c)
        try:
            rad.new_entry (c, entry.keys(), entry.values())
            conn.commit()
            return 0
        except sqlite3.IntegrityError:
            print ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
            return -1

def insert_product_entry (conn, product_name, init_product_table, entry):
    with conn:
        c = conn.cursor()
        p = init_product_table (product_name)
        p.create(c)
        keys = p.field_defs.keys()
        key_values = [entry[k] for k in keys]
        try:
            p.new_entry (c, keys, key_values)
            conn.commit()
            return 0
        except sqlite3.IntegrityError:
            print ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
            return -1

def insert_radiance_product_entry (conn, product_name, entry):
    return insert_product_entry (conn, product_name, init_radiance_product_table, entry)

def insert_other_product_entry (conn, product_name, entry):
    return insert_product_entry (conn, product_name, init_other_product_table, entry)

def get_radiance_keys (nc, keys):
    mirror_step = nc.variables["mirror_step"][:]
    keys["mirror_pos_beg"] = mirror_step.min()
    keys["mirror_pos_end"] = mirror_step.max()
    attr = nc.__dict__
    for k in Radiance_File_Attributes:
        if k in attr:
            keys[k] = attr[k]
    return keys

def remove_dot_prefix (name):
    if name.startswith('.'):
        return name[1:]
    else:
        return name

def process_file (conn, filename):

    basename = os.path.basename (filename)
    tok = basename.split('_')

    dont_process = ["tie"]
    if (tok[1] in dont_process):
        return 0

    product_name = '{}_{}'.format(tok[1], tok[2])

    nc = NetCDFFile(filename, "r")
    attr = nc.__dict__

    if not "time_coverage_start_since_epoch" in attr:
        print ("WARNING: missing attribute time_coverage_start_since_epoch; file={}".format (filename))
        return -1

    if (product_name == 'RAD_L1'):
        if attr["inr_status"] == "2":
            product_name = product_name + 'b'
        else:
            product_name = product_name + 'a'

    keys = {}
    keys["istart"] = int(attr["time_coverage_start_since_epoch"])
    keys["filename"] = remove_dot_prefix (basename)
    st = os.stat(filename)
    keys["mtime"] = st.st_mtime
    keys["size"] = st.st_size

    if (product_name in Radiance_Files):
        get_radiance_keys (nc, keys)
        status = insert_radiance_entry (conn, product_name, keys)
    elif (product_name in Radiance_Derived_Files):
        status = insert_radiance_product_entry (conn, product_name, keys)
    else:
        status = insert_other_product_entry (conn, product_name, keys)

    if status < 0:
        print('ERROR: processing file {}'.format(filename))

    return status

def make_db_path (arch_dir):
    db_basename = date.today().strftime("production_%Y%m.db")
    db_dir = os.path.join (arch_dir, "registry")
    if not os.path.isdir(db_dir):
        os.makedirs(db_dir)
    return os.path.join (db_dir, db_basename)

def register_files (db_path, filenames):

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")

    # Operate only on symbolic links!

    for fn in filenames:
        if (os.path.islink(fn)):
            status = process_file (conn, fn)
            if (status != 0):
                print('Error processing file: {}'.format(fn))
            os.remove(fn)

    conn.close()

def main():

    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if (arch_dir == None):
        printf ('*** Error: SDPC_ARCHIVE_DIR is not set')
        sys.exit(1)

    dir = os.path.join (arch_dir, 'registry/incoming')
    if not os.path.isdir (dir):
        os.makedirs(dir)

    db_path = make_db_path (arch_dir)

    while True:
        filenames = [os.path.join(dir,f) for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
        register_files (db_path, filenames)
        time.sleep (30)

    #sys.exit(status)

if __name__ == "__main__":
    main()
