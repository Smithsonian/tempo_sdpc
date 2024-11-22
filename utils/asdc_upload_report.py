#! /usr/bin/env python3

import os, sys
import time
from datetime import date
import dateutil.parser
import sqlite3
import re
import argparse

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3, "defer":100}
Asdc_Status_Lookup = {value: key for key, value in Asdc_Status.items()}

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def connect_database (db_file_path):
    conn = sqlite3.connect ("file:{}?mode=ro".format(db_file_path), uri=True)
    conn.execute("pragma foreign_keys=on")
    #conn.set_trace_callback(print)
    return conn

def get_product_table_names (cur):
    cur.execute ("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'");
    table_names = [item for t in cur.fetchall() for item in t]
    if table_names is None:
        return None
    # 'RAW' doesn't fit the table name pattern.
    append_raw = 'RAW' in table_names
    if append_raw:
        table_names.remove('RAW')
    table_names.sort(key = lambda x: x.split("_")[1])
    if append_raw:
        table_names.append('RAW')
    return table_names

def table_stats (cur, table_name, tbeg, tend):
    stats = {}
    # mean age at upload for files uploaded from this table in the specified time range
    sql = "select avg(asdc_upload_time-mtime) from {table_name} where asdc_upload_time > {tbeg} and asdc_upload_time < {tend}".format(**locals())
    cur.execute (sql)
    mean_upload_time = cur.fetchone()[0]
    if mean_upload_time is None:
        mean_upload_time = 0.0
    stats["mean_upload_time_min"] = mean_upload_time / 60.0

    # mean age at ingest for files ingested from this table in the specified time range
    sql = "select avg(asdc_ingest_time-mtime) from {table_name} where asdc_ingest_time > {tbeg} and asdc_ingest_time < {tend}".format(**locals())
    cur.execute (sql)
    mean_ingest_time = cur.fetchone()[0]
    if mean_ingest_time is None:
        mean_ingest_time = 0.0
    stats["mean_ingest_time_min"] = mean_ingest_time / 60.0

    # mean, and total size of files added to this table within in the specified time range
    sql = "select avg(size),sum(size) from {table_name} where mtime > {tbeg} and mtime < {tend}".format(**locals())
    cur.execute (sql)
    results = [item if item is not None else 0.0 for t in cur.fetchall() for item in t]
    stats["mean_size_MB"] = results[0] / 1.e6
    stats["total_size_GB"] = results[1] / 1.e9

    return stats

def table_status_summary (cur, table_name, tbeg, tend):
    sql = "select count(*) from {table_name} where mtime > {tbeg} and mtime < {tend}".format(**locals())
    cur.execute (sql)
    result = cur.fetchone()
    if result is None:
        num_files = 0
    else:
        num_files = result[0]
    sql = "select asdc_status,count(*) from {table_name} where mtime > {tbeg} and mtime < {tend} group by asdc_status".format(**locals())
    cur.execute (sql)
    results = cur.fetchall()
    # status_count = dict with Asdc_Status keys, with values initialized to zero
    status_count = dict.fromkeys(Asdc_Status, 0)
    for r in results:
        if r[0] in Asdc_Status_Lookup:
            status_count[Asdc_Status_Lookup[r[0]]] = r[1]
        else:
            eprint('*** Error: invalid asdc_status value: {}'.format(r[0]))
            sys.exit(1)
    return status_count, num_files

def print_table_summaries (db_file_path, table_list, tbeg, tend):
    print ("#")
    print (f"# pipeline: {os.environ['SDPC_PIPE_NAME']}")
    print (f"#   dbfile: {db_file_path}")
    if tbeg > 0:
        print("# Start time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tbeg))))
    print("#   End time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tend))))
    print ("#")
    with connect_database(db_file_path) as conn:
        cur = conn.cursor()
        if table_list is None:
            table_list = get_product_table_names (cur)
        print ("#       Product        mean        mean       mean    total    files     files   ingest   ingest")
        print ("#         table  age_upload  age_ingest       size     size    total  accepted  pending  problem")
        print ("#                     [min]       [min]       [MB]     [GB]                                     ")
        for table in table_list:
            status_count, num_files = table_status_summary (cur, table, tbeg, tend)
            stats = table_stats (cur, table, tbeg, tend)
            print ("%15s  %10.1f  %10.1f %9.1f %9.1f  %7d   %7d  %7d  %7d" % (table,
                   stats["mean_upload_time_min"],
                   stats["mean_ingest_time_min"],
                   stats["mean_size_MB"],
                   stats["total_size_GB"],
                   num_files,
                   status_count["accepted"],
                   status_count["uploaded"],
                   status_count["problem"]))

def main():
    parser = argparse.ArgumentParser(description='Summarize ASDC uploads')
    parser.add_argument('--dbfile', metavar='DBFILE', default=None,
                        help="sqlite database path")
    parser.add_argument('--start', default=None,
                        help="Start time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
    parser.add_argument('--end', default=None,
                        help="End time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
    parser.add_argument('--tables', metavar='TABLE', default=None, nargs="*",
                        help="Table name selection")
    # if len(sys.argv)==1:
    #     parser.print_usage(sys.stderr)
    #     sys.exit(0)
    args = parser.parse_args()

    if args.start is None:
        tbeg = 0
    else:
        tbeg = int(dateutil.parser.isoparse(args.start).timestamp())

    if args.end is None:
        tend = int(time.time())
    else:
        tend = int(dateutil.parser.isoparse(args.end).timestamp())

    if args.dbfile is None:
        db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
        if db_file_path == None:
            eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
            sys.exit(1)
    else:
        db_file_path = args.dbfile

    if not os.path.isfile (db_file_path):
        eprint ('nonexistent database file: {}'.format(db_file_path))
        sys.exit(0)

    print_table_summaries (db_file_path, args.tables, tbeg, tend)

if __name__ == "__main__":
    main()
