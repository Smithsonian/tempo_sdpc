#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import signal
import time
from threading import Event

import sqlite3
from subprocess import check_output
from netCDF4 import Dataset as NetCDFFile
import dateutil.parser as dp

Radiance_Files = ["RAD_L1a", "RAD_L1"]
Radiance_Products = ["CLDRR", "CLDO4", "BRO", "CHOCHO", "HCHO", "H2O", "NO2", "O3TOT", "O3PROF"]
Radiance_Derived_Files = [s + "_L2" for s in Radiance_Products] \
                       + [s + "_L3" for s in Radiance_Products]

Coverage_Time_Attributes = ["time_coverage_start_since_epoch", "time_coverage_end_since_epoch"]
Radiance_File_Attributes = Coverage_Time_Attributes \
                         + ["scan_num", "scan_type", "granule_num"]

Prefix = "register:"

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stderr, **kwargs)

def logprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stdout, **kwargs)

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
    fields["istart"] = "integer not null"
    fields["time_coverage_start_since_epoch"] = "float not null"
    fields["time_coverage_end_since_epoch"] = "float not null"
    fields["filename"] = "text"
    fields["path"] = "text"
    fields["mtime"] = "integer"
    fields["size"] = "integer"
    fields["versionid"] = "integer"
    fields["trend_status"] = "integer not null"
    fields["asdc_status"] = "integer"
    fields["asdc_status_met"] = "integer"
    fields["asdc_time_accepted"] = "integer"

def init_radiance_table (table_name):
    fields = {}
    define_common_fields (fields)
    fields["scan_type"] = "integer not null"
    fields["scan_num"] = "integer not null"
    fields["scan_id"] = "integer not null"
    fields["granule_num"] = "integer not null"
    fields["mirror_pos_beg"] = "integer not null"
    fields["mirror_pos_end"] = "integer not null"
    quals = "primary key(istart), unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_radiance_product_table (table_name):
    fields = {}
    define_common_fields (fields)
    fields["scan_id"] = "integer not null"
    quals = "unique(istart), foreign key (istart) references {}(istart)".format('RAD_L1')
    return Table_Type(table_name, fields, quals)

def init_dark_product_table (table_name):
    fields = {}
    define_common_fields (fields)
    fields["mean_exposure_time_per_coadd"] = "float"
    fields["mean_fpa_temp"] = "float"
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_other_product_table (table_name):
    fields = {}
    define_common_fields (fields)
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_raw_file_table (table_name):
    fields = {}
    fields["filename"] = "text"
    fields["path"] = "text"
    fields["istart"] = "integer not null"
    fields["mtime"] = "integer"
    fields["size"] = "integer"
    fields["asdc_status"] = "integer"
    fields["asdc_time_accepted"] = "integer"
    fields["asdc_status_met"] = "integer"
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def insert_raw_entry (conn, table_name, entry):
    c = conn.cursor()
    raw = init_raw_file_table(table_name)
    raw.create(c)
    try:
        raw.new_entry (c, entry.keys(), entry.values())
        return 0
    except sqlite3.IntegrityError:
        eprint ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
        return -1

def insert_radiance_entry (conn, table_name, entry):
    c = conn.cursor()
    rad = init_radiance_table(table_name)
    rad.create(c)
    try:
        rad.new_entry (c, entry.keys(), entry.values())
        conn.commit()
        return 0
    except sqlite3.IntegrityError:
        eprint ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
        return -1

def insert_product_entry (conn, product_name, init_product_table, entry):
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
        eprint ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
        return -1

def insert_radiance_product_entry (conn, product_name, entry):
    status = insert_product_entry (conn, product_name, init_radiance_product_table, entry)
    if status == 0:
        maybe_handle_scan_completion (conn, product_name, entry["scan_id"])
    return status

def insert_dark_product_entry (conn, product_name, entry):
    status = insert_product_entry (conn, product_name, init_dark_product_table, entry)
    return status

def insert_other_product_entry (conn, product_name, entry):
    status = insert_product_entry (conn, product_name, init_other_product_table, entry)
    return status

def get_dark_keys (nc, keys):
    variables_to_average = ["exposure_time_per_coadd", "fpa_temp"]
    for vname in variables_to_average:
        v = nc.variables[vname][:]
        keys["mean_{}".format(vname)] = v.mean()
    attr = nc.__dict__
    for k in Coverage_Time_Attributes:
        if k in attr:
            keys[k] = attr[k]
    return keys

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

def get_scan_id (path):
    scan_id = check_output (["level1_info", "-s", path])
    return int(scan_id)

def make_l3_path (l2_path):
    # l3_scan_dir
    product_dir = os.path.dirname(l2_path)
    granule_dir = os.path.dirname(product_dir)
    scan_dir = os.path.dirname (granule_dir)
    l3_scan_dir = scan_dir.replace ('/L2/', '/L3/')
    # l3_filename
    l3_filename = os.path.basename(l2_path)
    g_index = l3_filename.rfind('G')
    l3_filename = l3_filename[0:g_index].replace ('_L2_', '_L3_')
    # l3_path
    l3_path = os.path.join (l3_scan_dir, l3_filename) + '.nc'
    return l3_path

def write_pathlist (file_basename, pathlist):
    pipe_dir = os.getenv ("SDPC_RUN_DIR_MASTER")
    target_dir = os.path.join (pipe_dir, "stage/scans")
    os.makedirs (target_dir, exist_ok=True)
    # Ensure that the output file appears atomicly
    # write file to hidden_path, then rename to target_path
    hidden_path = os.path.join (target_dir, '.' + file_basename)
    target_path = os.path.join (target_dir, file_basename)
    fp = open (hidden_path, 'w')
    fp.write (pathlist)
    fp.close ()
    os.rename (hidden_path, target_path)

def collect_granule_paths (cur, product_name, scan_id):
    cur.execute ("select path from {} where scan_id == {} order by path".format(product_name, scan_id));
    paths = [item for t in cur.fetchall() for item in t]
    return paths

def count_archived_granules (cur, product_name, scan_id):
    cur.execute ("select count(scan_id) from {} where scan_id == {}".format(product_name, scan_id));
    return cur.fetchone()[0]

def handle_complete_scan (cur, product_name, scan_id):
    paths = collect_granule_paths (cur, product_name, scan_id)
    l3_path = make_l3_path (paths[0])
    pathlist = """product_name={}
l3_path={}
read -r -d '' l2_paths <<'EOF'
{}
EOF
""".format (product_name, l3_path, '\n'.join(paths))
    write_pathlist (os.path.basename (paths[0]), pathlist)

def maybe_handle_scan_completion (conn, product_name, scan_id):
    if not "L2" in product_name:
        return
    # This approach assumes that all of the RAD_L1a files will have been archived
    # by the time the associated L2 products start appearing.  This should be
    # true for normal processing because no L2 product generation starts until the
    # INR 2nd pass is finished, and the INR 2nd pass cannot finish until after all
    # of the RAD_L1a files have been delivered to INR.
    # However, this condition may not hold if near-real-time processing occurs
    # without the INR 2nd pass.
    cur = conn.cursor()
    num_radiance = count_archived_granules (cur, "RAD_L1a", scan_id)
    num_products = count_archived_granules (cur, product_name, scan_id)
    if num_products == num_radiance:
        handle_complete_scan (cur, product_name, scan_id)

def process_file (conn, filename, nc):

    basename = os.path.basename (filename)
    tok = basename.split('_')
    product_name = '{}_{}'.format(tok[1], tok[2])
    versionid = int(tok[3].strip('V'))

    attr = nc.__dict__

    if not "time_coverage_start_since_epoch" in attr:
        eprint ("WARNING: missing attribute time_coverage_start_since_epoch; file={}".format (filename))
        return -1

    # We use a count of the RAD_L1a granules to determine the number
    # of each Level 2 product type to expect from each scan.
    # That number of Level 2 products then triggers end-of-scan processing
    # for each Level 2 product type e.g. by L2_split and L2_regrid
    if product_name == 'RAD_L1' and attr["inr_status"] != "2":
        product_name = product_name + 'a'

    final_basename = remove_dot_prefix (basename)
    final_path = os.readlink (filename)
    dirname = os.path.dirname (final_path)
    st = os.stat(filename)

    # define common keys
    keys = {}
    keys["filename"] = final_basename
    keys["path"]     = final_path
    keys["size"]     = st.st_size
    keys["mtime"]    = int(st.st_mtime)
    keys["istart"]   = int(attr["time_coverage_start_since_epoch"])
    keys["versionid"] = versionid
    keys["trend_status"] = 0
    keys["asdc_status"] = 0
    keys["asdc_time_accepted"] = 0
    keys["time_coverage_start_since_epoch"] = attr["time_coverage_start_since_epoch"]
    keys["time_coverage_end_since_epoch"] = attr["time_coverage_end_since_epoch"]

    # Look for a .met file
    if os.path.exists(final_path + ".met"):
        keys["asdc_status_met"] = 0   # new
    else:
        keys["asdc_status_met"] = -2  # nonexistent

    if product_name in Radiance_Files:
        keys["scan_id"] = get_scan_id (final_path)
        get_radiance_keys (nc, keys)
        status = insert_radiance_entry (conn, product_name, keys)
    elif product_name in Radiance_Derived_Files:
        keys["scan_id"] = get_scan_id (final_path)
        status = insert_radiance_product_entry (conn, product_name, keys)
    elif product_name == "DRK_L1":
        get_dark_keys (nc, keys)
        status = insert_dark_product_entry (conn, product_name, keys)
    else:
        status = insert_other_product_entry (conn, product_name, keys)

    if status < 0:
        eprint('ERROR: processing file {}'.format(filename))
    else:
        logprint ('{}: {}'.format(product_name, final_basename), flush=True)

    return status

def process_file_raw (conn, filename):

    basename = os.path.basename (filename)
    final_path = os.readlink (filename)
    st = os.stat (filename)

    # TEMPO_GRDDP_YYYYMMDDThhmmssZ_Cccccc.raw
    # TEMPO_CCSDS_YYYYMMDDThhmmssZ_Cccccc.raw
    tok = basename.split('_')
    table_name = tok[1];
    istart = int(dp.parse(tok[2]).timestamp())

    keys = {}
    keys["filename"] = basename
    keys["path"] = final_path
    keys["istart"] = istart
    keys["mtime"] = int(st.st_mtime)
    keys["size"] = st.st_size
    keys["asdc_status"] = 0
    keys["asdc_time_accepted"] = 0
    keys["asdc_status_met"] = -2 # nonexistent

    status = insert_raw_entry (conn, table_name, keys)

    if status < 0:
        eprint('ERROR: processing file {}'.format(filename))
    else:
        logprint ('{}: {}'.format(table_name, basename), flush=True)

    return status

def connect_database (db_path):
    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")
    return conn

def register_files (db_path, filenames):

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    # Operate only on symbolic links!

    for fn in filenames:
        if os.path.islink(fn):
            with connect_database (db_path) as conn:
                if fn.endswith ('.nc'):
                    with NetCDFFile (fn, "r") as nc:
                        status = process_file (conn, fn, nc)
                elif fn.endswith ('.raw'):
                    status = process_file_raw (conn, fn)
                else:
                    status = -1
            if status != 0:
                eprint('Error processing file: {}'.format(fn))
            os.remove(fn)

def collect_filenames (dir):
    filenames = []
    for f in os.listdir(dir):
        path = os.path.join(dir, f)
        if os.path.isfile(path):
            filenames.append (path)
    return filenames

class Registry:
    def __init__ (self, incoming_dir, file_path):
        self.incoming_dir = incoming_dir
        self.file_path = file_path

def init_registry ():
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if arch_dir == None:
        eprint ('*** Error: SDPC_ARCHIVE_DIR is not set')
        sys.exit(1)

    incoming_dir = os.path.join (arch_dir, 'registry/incoming')
    if not os.path.isdir (incoming_dir):
        os.makedirs(incoming_dir)

    return Registry (incoming_dir, db_file_path)

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGHUP, self.handler)
    signal.signal(signal.SIGTERM, self.handler)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

def main():

    reg = init_registry()
    sig = Signal_Catcher()

    logprint ("Started", flush=True)

    while not sig.caught():
        filenames = collect_filenames (reg.incoming_dir)
        if len(filenames) > 0:
            register_files (reg.file_path, filenames)
        sig.wait(10)

    logprint ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
