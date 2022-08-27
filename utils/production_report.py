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

def get_matching_table_names (cur, pat):
    cur.execute ("select name from sqlite_master where type = 'table' and name not like 'sqlite_%'");
    table_names = [item for t in cur.fetchall() for item in t]
    if table_names is None:
        return None
    return [item for item in table_names if pat in item]

def level1_rad_stats (cur, tbeg, tend):
    rad_tables = get_matching_table_names (cur, 'RAD')
    if rad_tables is None:
        return None
    if 'RAD_L0' not in rad_tables:
        return None
    if 'RAD_L1' not in rad_tables:
        return None
    stats = {}
    # time elapsed between creation of RAD_L0 and completion of RAD_L1
    sql = 'select avg(RAD_L1.mtime-RAD_L0.mtime) from '\
          'RAD_L1 inner join RAD_L0 on RAD_L0.istart = RAD_L1.istart' \
          'where RAD_L0.mtime > {tbeg} and RAD_L0.mtime < {tend}'.format(**locals())
    cur.execute (sql)
    stats["mean_elapsed_sec"] = cur.fetchone()[0]
    return stats

def level2_stats (cur, table_name, tbeg, tend):
    # time elapsed between completion of L1 radiance file and completion of specified L2 product
    stats = {}
    sql = 'select avg({table_name}.mtime-RAD_L1.mtime) from '\
          '{table_name} inner join RAD_L1 on RAD_L1.istart = {table_name}.istart ' \
          'where RAD_L1.mtime > {tbeg} and RAD_L1.mtime < {tend}'.format(**locals())
    cur.execute (sql)
    mean_elapsed_sec = cur.fetchone()[0]
    if mean_elapsed_sec is None:
        mean_elapsed_sec = 0.0
    stats["mean_elapsed_sec"] = mean_elapsed_sec
    return stats

def print_latencies (tbeg, tend):
    print ("#")
    print (f"#   SDPC_PIPE_NAME: {os.environ['SDPC_PIPE_NAME']}")
    if tbeg > 0:
        print("# Start time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tbeg))))
    print("#   End time: %s" % (time.strftime ('%Y-%m-%dT%H:%M:%SZ', time.gmtime(tend))))
    print ("#")
    with connect_database() as conn:
        cur = conn.cursor()
        print ("#       Product        mean")
        print ("#         table  production")
        print ("#                     [min]")
        stats = level1_rad_stats (cur, tbeg, tend)
        if stats is not None:
            print ("%15s  %10.1f" % ('RAD_L1', stats["mean_elapsed_sec"]/60.0))
        level2_tables = get_matching_table_names (cur, "L2")
        for table in level2_tables:
            stats = level2_stats (cur, table, tbeg, tend)
            print ("%15s  %10.1f" % (table, stats["mean_elapsed_sec"]/60.0))

def main():
    parser = argparse.ArgumentParser(description='Summarize production latencies')
    parser.add_argument('--start', default=None,
                        help="Start time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
    parser.add_argument('--end', default=None,
                        help="End time, ISO format e.g. YYYY-MM-DDThh:mm:ss[Z]")
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

    print_latencies (tbeg, tend)

if __name__ == "__main__":
    main()
