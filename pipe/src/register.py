#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import traceback

import os, sys
import signal
import time
from threading import Event

import sqlite3
import argparse
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

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3, "defer":100}

Prefix = "register:"

Have_Rad_L1_Table = False
Have_Rad_L1a_Table = False
Alt_Rad_L1a_Dbfile_Path = None
HCHO_Needs_Destripe = False

Reload_Config = False

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

def basic_fields_dict ():
    fields = {}
    fields["filename"] = "text"
    fields["path"] = "text"
    fields["istart"] = "integer not null"
    fields["mtime"] = "integer"
    fields["size"] = "integer"
    fields["asdc_status"] = "integer"
    fields["asdc_status_met"] = "integer"
    fields["asdc_upload_time"] = "integer"
    fields["asdc_ingest_time"] = "integer"
    fields["asdc_disposition"] = "text"
    return fields

def common_fields_dict():
    fields = basic_fields_dict()
    fields["time_coverage_start_since_epoch"] = "float not null"
    fields["time_coverage_end_since_epoch"] = "float not null"
    fields["versionid"] = "integer"
    fields["trend_status"] = "integer not null"
    return fields

def init_radiance_table (table_name):
    fields = common_fields_dict()
    fields["scan_type"] = "integer not null"
    fields["scan_num"] = "integer not null"
    fields["scan_id"] = "integer not null"
    fields["granule_num"] = "integer not null"
    fields["mirror_pos_beg"] = "integer not null"
    fields["mirror_pos_end"] = "integer not null"
    quals = "primary key(istart), unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_radiance_product_table (table_name):
    fields = common_fields_dict()
    fields["scan_id"] = "integer not null"
    if Have_Rad_L1_Table:
        quals = "unique(istart), foreign key (istart) references RAD_L1(istart)"
    else:
        quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_dark_product_table (table_name):
    fields = common_fields_dict()
    fields["mean_exposure_time_per_coadd"] = "float"
    fields["mean_fpa_temp"] = "float"
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_other_product_table (table_name):
    fields = common_fields_dict()
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_raw_file_table (table_name):
    fields = basic_fields_dict()
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def init_corrfile_table (table_name):
    fields = basic_fields_dict()
    fields["tstart"] = "integer not null"
    fields["tend"] = "integer not null"
    fields["begin_hour_utc"] = "float"
    fields["end_hour_utc"] = "float"
    fields["num_mirror_pos"] = "integer"
    quals = "unique(istart)"
    return Table_Type(table_name, fields, quals)

def insert_raw_entry (conn, table_name, entry):
    c = conn.cursor()
    raw = init_raw_file_table(table_name)
    raw.create(c)
    try:
        raw.new_entry (c, entry.keys(), entry.values())
        conn.commit()
        return 0
    except sqlite3.IntegrityError:
        eprint ('ERROR: duplicate primary key: istart={}'.format(entry["istart"]))
        return -1

def insert_corrfile_entry (conn, table_name, entry):
    c = conn.cursor()
    tbl = init_corrfile_table(table_name)
    tbl.create(c)
    try:
        tbl.new_entry (c, entry.keys(), entry.values())
        conn.commit()
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
    pipe_dir = os.getenv ("SDPC_PIPE_DIR")
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
    """
    To see when all products for a given scan are complete, we compare the number
    of products with the number of archived RAD_L1a files for that scan.

    During normal processing, all of the RAD_L1a files will have been archived
    by the time the associated L2 products start appearing because no L2 product
    generation starts until the INR 2nd pass is finished, and the INR 2nd pass
    cannot finish until after all of the RAD_L1a files have been delivered to INR.

    During L2 reprocessing, the RAD_L1a files may be in a different archive, so we
    may need to query a different sqlite database for a count of RAD_L1a files.

    If we support near-real-time processing, then the INR 2nd pass may not be performed.
    In that case, it may be possible for L2 products to be produced before all RAD_L1a
    products from a given scan are archived.  In that case, we may need a different
    method to detect scan completion for a given NRT product.
    """
    if not "L2" in product_name:
        return
    cur = conn.cursor()
    num_products = count_archived_granules (cur, product_name, scan_id)

    if Have_Rad_L1a_Table:
        num_radiance = count_archived_granules (cur, "RAD_L1a", scan_id)
    elif Alt_Rad_L1a_Dbfile_Path == None:
        return
    else:
        with connect_database (Alt_Rad_L1a_Dbfile_Path) as conn_sep:
            num_radiance = count_archived_granules (conn_sep.cursor(), "RAD_L1a", scan_id)

    if num_products == num_radiance:
        handle_complete_scan (cur, product_name, scan_id)

def table_exists (conn, table_name):
    cur = conn.cursor()
    cur.execute ("SELECT name FROM sqlite_master WHERE type='table' AND name='{}';".format(table_name))
    result = cur.fetchone()
    return result != None

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
    keys["asdc_status"] = Asdc_Status["new"]
    keys["asdc_upload_time"] = 0
    keys["asdc_ingest_time"] = 0
    keys["asdc_disposition"] = ""
    keys["time_coverage_start_since_epoch"] = attr["time_coverage_start_since_epoch"]
    keys["time_coverage_end_since_epoch"] = attr["time_coverage_end_since_epoch"]

    # Look for a .met file
    if os.path.exists(final_path + ".met"):
        keys["asdc_status_met"] = Asdc_Status["new"]
    else:
        keys["asdc_status_met"] = Asdc_Status["nonexistent"]

    # Some products require additional processing steps before
    # uploading to ASDC.  Such products are initially registered
    # with status "defer", which is later updated to "new" (elsewhere)
    # upon completion of the final processing step.
    # currently: NO2_L2 waits for strat/trop separation
    #            HCHO_L2 may wait for destriping/background correction
    if product_name == 'NO2_L2':
        defer_asdc_upload = True
    elif (product_name == 'HCHO_L2' and HCHO_Needs_Destripe):
        defer_asdc_upload = (('destriping_correction' not in nc['support_data'].variables) and
                             ('background_correction' not in nc['support_data'].variables))
    else:
        defer_asdc_upload = False

    if defer_asdc_upload:
        keys["asdc_status"] = Asdc_Status["defer"]
        keys["asdc_status_met"] = Asdc_Status["defer"]

    global Have_Rad_L1a_Table
    Have_Rad_L1a_Table = table_exists (conn, 'RAD_L1a')
    global Have_Rad_L1_Table
    Have_Rad_L1_Table = table_exists (conn, 'RAD_L1')

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

    # tempo_YYYYMMDDThhmmssZ_<suffix>.tar
    tok = basename.split('_')
    table_name = "RAW";
    istart = int(dp.parse(tok[1]).timestamp())

    keys = {}
    keys["filename"] = basename
    keys["path"] = final_path
    keys["istart"] = istart
    keys["mtime"] = int(st.st_mtime)
    keys["size"] = st.st_size
    keys["asdc_status"] = Asdc_Status["new"]
    keys["asdc_upload_time"] = 0
    keys["asdc_ingest_time"] = 0
    keys["asdc_status_met"] = Asdc_Status["nonexistent"]
    keys["asdc_disposition"] = ""

    status = insert_raw_entry (conn, table_name, keys)

    if status < 0:
        eprint('ERROR: processing file {}'.format(filename))
    else:
        logprint ('{}: {}'.format(table_name, basename), flush=True)

    return status

def convert_hms_to_float (hms_array):
    hour_f = []
    for hms in hms_array:
        fields = hms.split(':')
        h_f = int(fields[0]) + (int(fields[1]) + int(fields[2])/60.0)/60.0
        hour_f.append (h_f)
    return hour_f

def process_file_corr (conn, filename, nc):

    basename = os.path.basename (filename)
    final_path = os.readlink (filename)
    st = os.stat (filename)

    # Example:   TEMPO_RADREF_L1_V01_YYYYMMDD_S123456789_E123456789_S001.nc
    # Example: TEMPO_DSTRHCHO_L2_V01_YYYYMMDD_S123456789_E123456789_S001.nc
    tok = basename.split('_')
    table_name = "_".join ([tok[1], tok[2]])
    tstart = int(tok[5].strip('S'))
    tend = int(tok[6].strip('E'))

    keys = {}
    keys["filename"] = basename
    keys["path"] = final_path
    keys["istart"] = tstart
    keys["mtime"] = int(st.st_mtime)
    keys["size"] = st.st_size
    keys["asdc_status"] = Asdc_Status["new"]
    keys["asdc_upload_time"] = 0
    keys["asdc_ingest_time"] = 0
    keys["asdc_status_met"] = Asdc_Status["nonexistent"]
    keys["asdc_disposition"] = ""

    begin_time = nc.getncattr('begin_time')
    end_time = nc.getncattr('end_time')
    begin_time = convert_hms_to_float (begin_time.split(','))
    end_time = convert_hms_to_float (end_time.split(','))

    keys["tstart"] = tstart
    keys["tend"] = tend
    keys["begin_hour_utc"] = min(begin_time)
    keys["end_hour_utc"] = max(end_time)
    keys["num_mirror_pos"] = nc.getncattr("num_mirror_pos")

    status = insert_corrfile_entry (conn, table_name, keys)

    if status < 0:
        eprint('ERROR: processing file {}'.format(filename))
    else:
        logprint ('{}: {}'.format(table_name, basename), flush=True)

    return status

def connect_database (db_path):
    conn = sqlite3.connect (db_path)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def register_one_file (db_path, fn):
    status = -1
    with connect_database (db_path) as conn:
        basename = os.path.basename(fn)
        if basename.startswith ('TEMPO_RADREF') or basename.startswith('TEMPO_DSTR'):
            with NetCDFFile (fn, "r") as nc:
                status = process_file_corr (conn, fn, nc)
        elif fn.endswith ('.nc'):
            with NetCDFFile (fn, "r") as nc:
                status = process_file (conn, fn, nc)
        elif fn.endswith ('.tar'):
            status = process_file_raw (conn, fn)
        if status != 0:
            eprint('Error processing file: {}'.format(fn))
    return status

def move_failing (fn):
    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if arch_dir == None:
        eprint ('*** Error: SDPC_ARCHIVE_DIR is not set')
        return
    fail_dir = os.path.join (arch_dir, 'registry/failed')
    if not os.path.isdir (fail_dir):
        os.makedirs(fail_dir)
    new_path = os.path.join (fail_dir, os.path.basename(fn))
    logprint ('moving {} to {}'.format(fn, new_path))
    os.rename (fn, new_path)

def register_files (db_path, filenames):

    # For back-compatibility sqlite has foreign keys turned off by default,
    # and foreign_keys=off is ALWAYS stored in the database, regardless of
    # the runtime setting when the database was created.  For this reason,
    # we apparently need to turn it on explicitly, each time the database
    # connection is established.

    # Operate only on symbolic links!

    status_list = []

    for fn in filenames:
        if os.path.islink(fn):
            try:
                status = register_one_file (db_path, fn)
            except BaseException as e:
                print('An exception occurred: {}'.format(e))
                print(traceback.print_exc())
                status_list.append(-1)
                move_failing(fn)
            else:
                status_list.append(status)
                os.remove(fn)

    return status_list

def collect_filenames (dir):
    entries = (os.path.join (dir, f) for f in os.listdir(dir))
    entries = sorted(entries, key=os.path.getmtime)
    filenames = []
    for path in entries:
        if os.path.isfile(path):
            filenames.append (path)
    return filenames

class Registry:
    def __init__ (self, incoming_dir, file_path):
        self.incoming_dir = incoming_dir
        self.file_path = file_path

def check_hcho_destripe_config():
    # Does HCHO_L2 receive a destriping/background correction?
    global HCHO_Needs_Destripe
    setting = check_output (["config_setting", "destripe.HCHO.apply"])
    HCHO_Needs_Destripe = setting != 0

def load_config():
    global Reload_Config
    Reload_Config = False
    check_hcho_destripe_config()

def init_registry ():
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    db_l1_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE_L1")
    if db_l1_file_path != None and db_l1_file_path != db_file_path:
        global Alt_Rad_L1a_Dbfile_Path
        Alt_Rad_L1a_Dbfile_Path = db_l1_file_path

    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if arch_dir == None:
        eprint ('*** Error: SDPC_ARCHIVE_DIR is not set')
        sys.exit(1)

    incoming_dir = os.path.join (arch_dir, 'registry/incoming')
    if not os.path.isdir (incoming_dir):
        os.makedirs(incoming_dir)

    load_config()

    return Registry (incoming_dir, db_file_path)

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGTERM, self.handler)
    signal.signal(signal.SIGHUP, self.config_update)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

  def config_update(self,signum,frame):
      global Reload_Config
      Reload_Config = True
      self.signum = signum

def run_as_service (reg):

    sig = Signal_Catcher()

    logprint ("Started", flush=True)

    while not sig.caught():
        if Reload_Config:
            logprint ('Caught signal = SIGHUP: updating configuration', flush=True)
            load_config()
        filenames = collect_filenames (reg.incoming_dir)
        if len(filenames) > 0:
            status_list = register_files (reg.file_path, filenames)
        sig.wait(10)

    logprint ("Exiting: caught signal = {}".format(sig.signum))

def main():
    parser = argparse.ArgumentParser(description='register data products in a sqlite database')
    parser.add_argument('--service', help="run as a service", action='store_true')
    parser.add_argument('filename', help="netCDF4 file name", nargs='*', default=None)
    args = parser.parse_args()

    reg = init_registry()

    if args.service:
        run_as_service (reg)
    elif args.filename != None:
        status_list = register_files (reg.file_path, args.filename)
        if any(status_list):
            sys.exit(1)

if __name__ == "__main__":
    main()
