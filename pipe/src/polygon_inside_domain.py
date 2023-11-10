#! /usr/bin/env python3

import sys
import os
import argparse
import csv

from netCDF4 import Dataset as NetCDFFile
import numpy as np

Verbose = False

def read_bpoly (filename):
    # parse OGS boundary polygon attribute
    with NetCDFFile(filename) as nc:
        geospatial_bounds = nc.getncattr("geospatial_bounds")
    bdry_pts = geospatial_bounds.lstrip('POLYGON((').rstrip('))').split(',')
    bdry_lon = [float(p.split(' ')[1]) for p in bdry_pts]
    bdry_lat = [float(p.split(' ')[0]) for p in bdry_pts]
    return bdry_lon, bdry_lat

def strip_comments(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

# strip comments and leading/trailing spaces
def read_csv_columns (path):
    with open (path,'r') as fp:
        reader = csv.reader (strip_comments(fp))
        headers = next (reader, None)
        headers = [h.strip() for h in headers]
        column = {}
        for h in headers:
            column[h] = []
        for row in reader:
            for h, v in zip(headers, row):
                column[h].append(v.strip())
    return column

# ray intersection method
def point_in_polygon (xx,yy,px,py):
    count = 0

    xa = px[0]
    ya = py[0]

    if xx == xa and yy == ya:
        return 1

    num = len(px)

    for i in range(1,num):
        xb = px[i]
        yb = py[i]

        if xx == xb and yy == yb:
            return 1

        # does line x=xx intersect this edge?
        if (xb - xx) * (xx - xa) < 0:
            xa=xb
            ya=yb
            continue

        if xa != xb:
            yi = ya + (xx - xa) * (yb - ya) / (xb - xa)
            # intersection lies below point?
            if yi <= yy:
                count += 1
        else:
            # point lies on vertical edge?
            if (yb - yy) * (yy - ya) >= 0:
                return 1
        xa=xb
        ya=yb

    # ray from an inside point intersects the polygon
    # in an odd number of points.
    if count % 2 != 0:
        return True
    else:
        return False

def polygon_inside_boundary (lons, lats, b_lons, b_lats):
    for p in zip (lons,lats):
        pt_inside = point_in_polygon (p[0], p[1], b_lons, b_lats)
        if pt_inside is False:
            if Verbose is True:
                print ("point is outside: {}".format(p))
            break
    return pt_inside

def main():
    parser = argparse.ArgumentParser(description='validate geospatial bounding polygon')
    parser.add_argument('--domain', help="CSV file defining domain polygon")
    parser.add_argument('--verbose', action="store_true")
    parser.add_argument('filename', help="data product filename")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    global Verbose
    Verbose = args.verbose

    if not os.path.isfile(args.domain):
        print ("*** Error: file not found: {}".format(args.domain))
        sys.exit(1)

    b = read_csv_columns (args.domain)
    b_lons = [float(x) for x in b["longitude"][:]]
    b_lats = [float(x) for x in b["latitude"][:]]

    if not os.path.isfile(args.filename):
        print ("*** Error: file not found: {}".format(args.filename))
        sys.exit(1)

    (lons, lats) = read_bpoly (args.filename)
    is_inside = polygon_inside_boundary (lons, lats, b_lons, b_lats)

    if is_inside is True:
        print ("yes")
    else:
        print ("no")

if __name__ == '__main__':
    main()

