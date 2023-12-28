#! /usr/bin/env python3

import os, sys
import glob
import time
from datetime import date
import dateutil.parser
import sqlite3
import re
import argparse

import numpy as np

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3, "defer":100}
Asdc_Status_Lookup = {value: key for key, value in Asdc_Status.items()}

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def choose_db_file_path():
    """
    Directly accessing SDPC_ARCHIVE_DBFILE competes with ongoing pipeline processes
    for database access and could cause a random failure due to database contention.
    Ideally, the pipeline should be robust enough so that this isn't a problem,
    (and we open the database file read-only, so surely that helps!)
    but until we're certain that this routine isn't causing contention, it's best to
    access the live database file only when it's essential to do so.  For this purpose
    of this script, querying a recent backup copy should be sufficient.
    """
    db_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    db_files = glob.glob (os.path.join (db_dir, "backup/archived_products.sqlite.*"))
    db_file_path = max (db_files, key=os.path.getmtime)
    if not os.path.isfile (db_file_path):
        eprint ('*** Error: finding sqlite database file path in {}'.format(db_dir))
        sys.exit(1)
    return db_file_path

def connect_database (dbfile):
    # QUESTION:  Does opening the database read-only avoid contention that might lead to a "database is locked" error?
    #conn = sqlite3.connect (dbfile)
    conn = sqlite3.connect ("file:{}?mode=ro".format(dbfile), uri=True)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def measurement_to_upload_latency (cur, table_name, tbeg, tend):
    tempo_epoch_timet = 315964800
    max_proc_delta = 86400
    accepted = Asdc_Status["accepted"]
    sql = 'select delta,cume_dist from ('\
                  'select (asdc_upload_time-(time_coverage_start_since_epoch+{tempo_epoch_timet})) as delta,'\
                                                                             'cume_dist() over win as cume_dist '\
                  'from {table_name} '\
                       'where asdc_status == {accepted} and mtime-(time_coverage_start_since_epoch+{tempo_epoch_timet}) < {max_proc_delta} and mtime between {tbeg} and {tend} '\
                  'window win as (order by asdc_upload_time-time_coverage_start_since_epoch)) t '.format(**locals())
                  # 'where cume_dist > {percentile} order by cume_dist limit 1'
    cur.execute (sql)
    rows = cur.fetchall()
    delta = []
    cume_dist = []
    for (a,b) in rows:
        delta.append(a)
        cume_dist.append(b)
    return {"delta":np.asarray(delta), "cume_dist":np.asarray(cume_dist), "num":len(delta)}

def main():
    parser = argparse.ArgumentParser(description='Summarize production latencies')
    parser.add_argument('--dbfile', default=None, help = "sqlite database file path")
    parser.add_argument('--start', default=None,
                        help="Start time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
    parser.add_argument('--end', default=None,
                        help="End time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
    # if len(sys.argv)==1:
    #     parser.print_usage(sys.stderr)
    #     sys.exit(0)
    args = parser.parse_args()

    # Establish time range
    if args.start is None:
        tbeg = 0
    else:
        tbeg = int(dateutil.parser.isoparse(args.start).timestamp())

    if args.end is None:
        tend = int(time.time())
    else:
        tend = int(dateutil.parser.isoparse(args.end).timestamp())

    # Get sqlite database file path
    dbfile = args.dbfile
    if dbfile is None:
        dbfile = choose_db_file_path()

    if not os.path.isfile (dbfile):
        eprint ("Invalid database path: {}".format(dbfile))
        sys.exit(1)

    print ("# dbfile = {}".format(dbfile))

    with connect_database (dbfile) as conn:
        cur = conn.cursor()
        cur.execute ("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'");
        table_names = [item for t in cur.fetchall() for item in t]

    excluded_tables = ["RAW", "RADREF_L1", "DSTRHCHO_L2"]

    product = []
    latency_50 = []
    latency_90 = []
    num = []
    for table in table_names:
        if "L0" in table or table in excluded_tables:
            continue
        with connect_database(dbfile) as conn:
            d = measurement_to_upload_latency (conn.cursor(), table, tbeg, tend)
            delta = d["delta"]
            if d["num"] < 2:
                continue
            product.append(table)
            num.append(d["num"])
            cume_dist = d["cume_dist"]
            t50 = np.interp (0.5, cume_dist, delta)
            t90 = np.interp (0.9, cume_dist, delta)
            latency_50.append(t50/3600.0)
            latency_90.append(t90/3600.0)

    indices = np.argsort(latency_50)
    print ("product_name,num,latency_50pcnt,latency_90pcnt")
    for i in indices:
        print("%s,%d,%f,%f" % (product[i], num[i], latency_50[i], latency_90[i]))

if __name__ == "__main__":
    main()
