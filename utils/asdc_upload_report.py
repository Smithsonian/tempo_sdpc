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

def connect_database ():
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    conn = sqlite3.connect (db_file_path)
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

def table_accepted_summary (cur, table_name, tbeg, tend):
    accepted = Asdc_Status["accepted"]
    sql = "select avg(asdc_upload_time-mtime),avg(asdc_ingest_time-mtime),avg(size),count(*),sum(size) from {table_name} where mtime > {tbeg} and mtime < {tend} and asdc_status == {accepted}".format(**locals())
    cur.execute (sql)
    results = [item for t in cur.fetchall() for item in t]
    return {'upload_time':results[0], 'ingest_time':results[1], 'size':results[2], 'count':results[3], 'total_size':results[4]}

def table_status_summary (cur, table_name, tbeg, tend):
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
    return status_count

def print_table_summaries (table_list, tbeg, tend):
    print ("#")
    if tbeg > 0:
        print("# Start time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tbeg))))
    print("#   End time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tend))))
    print ("#")
    with connect_database() as conn:
        cur = conn.cursor()
        if table_list is None:
            table_list = get_product_table_names (cur)
        print ("#  Product        mean        mean       mean    total     files   ingest   ingest")
        print ("#    table  age_upload  age_ingest       size     size  accepted  pending  problem")
        print ("#                [min]       [min]       [MB]     [GB]                            ")
        for table in table_list:
            status_count = table_status_summary (cur, table, tbeg, tend)
            summary = table_accepted_summary (cur, table, tbeg, tend)
            if summary["count"] == 0:
                summary = {'upload_time':-1, 'ingest_time':-1, 'size':-1, 'count':0, 'total_size':-1}
            print ("%10s  %10.1f  %10.1f %9.1f %9.1f   %7d  %7d  %7d" % (table,
                   summary["upload_time"]/60.0, summary["ingest_time"]/60.0, summary["size"]/1.e6,
                   summary["total_size"]/1.e9, summary["count"],
                   status_count["uploaded"], status_count["problem"]))

def tlimits_from_tbounds (tbounds):
    if tbounds is not None:
        tbeg_obj = dateutil.parser.isoparse(tbounds[0])
        tend_obj = dateutil.parser.isoparse(tbounds[1])
        tbeg = int(tbeg_obj.timestamp())
        tend = int(tend_obj.timestamp())
    else:
        tbeg = 0
        tend = int(time.time())
    return tbeg, tend

def main():
    parser = argparse.ArgumentParser(description='Summarize ASDC uploads')
    parser.add_argument('--interval', metavar=('BEGIN','END',), default=None, nargs=2,
                        help="Time interval selection")
    parser.add_argument('--tables', metavar='TABLE', default=None, nargs="*",
                        help="Table name selection")
    # if len(sys.argv)==1:
    #     parser.print_usage(sys.stderr)
    #     sys.exit(0)
    args = parser.parse_args()

    tbeg, tend = tlimits_from_tbounds (args.interval)

    print_table_summaries (args.tables, tbeg, tend)

if __name__ == "__main__":
    main()
