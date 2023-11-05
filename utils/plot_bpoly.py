#! /usr/bin/env python3

import sys
import glob
import os
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from netCDF4 import Dataset as NetCDFFile
import numpy as np

#+++ begin lines to stop shapely deprecation warnings
import shapely
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
#--- end lines to stop shapely deprecation warnings

import cartopy.crs as ccrs
import cartopy.feature

def read_bpoly (filename):
    # parse OGS boundary polygon attribute
    with NetCDFFile(filename) as nc:
        geospatial_bounds = nc.getncattr("geospatial_bounds")
    bdry_pts = geospatial_bounds.lstrip('POLYGON((').rstrip('))').split(',')
    bdry_lon = [float(p.split(' ')[1]) for p in bdry_pts]
    bdry_lat = [float(p.split(' ')[0]) for p in bdry_pts]
    return bdry_lon, bdry_lat

def main():
    parser = argparse.ArgumentParser(description='plot geospatial bounding polygon')
    parser.add_argument('--outfile', help="path to output file")
    parser.add_argument('files', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    fig = plt.figure()
    ax = plt.subplot (1,1,1, projection=ccrs.Miller())
    ax.set_extent ([-155, -30, 15, 62], ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5)
    gl.top_labels=False
    gl.right_labels=False
    gl.xlabel_style = {'size':6}
    gl.ylabel_style = {'size':6}
    ax.add_feature (cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature (cartopy.feature.BORDERS, linewidth=0.5)
    ax.add_feature (cartopy.feature.STATES, linewidth=0.5)

    for path in args.files:
        (lons, lats) = read_bpoly (path)
        ax.plot (lons, lats, lw=0.5, transform=ccrs.PlateCarree())

    #plt.show()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

