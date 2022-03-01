#! /usr/bin/env python3

import os, sys
import time
from datetime import date
import dateutil.parser
import sqlite3
import re
import argparse

Asdc_Status = {"nonexistent":-2, "problem":-1, "new": 0, "pending":1, "uploaded":2, "accepted":3}

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
    if table_names is not None:
        table_names.remove('GRDDP')
        table_names.sort(key = lambda x: x.split("_")[1])
    return table_names

def table_accepted_summary (cur, table_name, tbounds):
    accepted = Asdc_Status["accepted"]
    if tbounds is not None:
        tbeg_obj = dateutil.parser.isoparse(tbounds[0])
        tend_obj = dateutil.parser.isoparse(tbounds[1])
        tbeg = int(tbeg_obj.timestamp())
        tend = int(tend_obj.timestamp())
    else:
        tbeg = 0
        tend = int(time.time())
    sql = "select avg(asdc_upload_time-mtime),avg(asdc_ingest_time-mtime),avg(size),count(*) from {table_name} where mtime > {tbeg} and mtime < {tend} and asdc_status == {accepted}".format(**locals())

    cur.execute (sql)
    results = [item for t in cur.fetchall() for item in t]

    return {'upload_time':results[0], 'ingest_time':results[1], 'size':results[2], 'count':results[3]}

def print_table_summaries (table_list, tbounds):
    with connect_database() as conn:
        cur = conn.cursor()
        if table_list is None:
            table_list = get_product_table_names (cur)
        print ("#  Product        mean        mean       mean    file")
        print ("#    table  age_upload  age_ingest       size   count")
        print ("#                [min]       [min]       [MB]        ")
        for table in table_list:
            summary = table_accepted_summary (cur, table, tbounds)
            if summary["upload_time"] is None or summary["ingest_time"] is None:
                continue;
            print ("%10s  %10.1f  %10.1f %9.1f  %7d" % (table,
            summary["upload_time"]/60.0, summary["ingest_time"]/60.0, summary["size"]/1.e6,
            summary["count"]))

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

    print_table_summaries (args.tables, args.interval)

if __name__ == "__main__":
    main()
