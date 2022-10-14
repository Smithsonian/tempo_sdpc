#! /usr/bin/env python3

import sys
import os
import math

import matplotlib
matplotlib.use('Agg')

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

import numpy as np
from netCDF4 import Dataset as NetCDFFile
import argparse

def find_point (lons_deg, lats_deg, lon0_deg, lat0_deg):
    degtorad = math.pi/180.0
    lats = lats_deg * degtorad
    lons = lons_deg * degtorad
    ny,nx = lats.shape
    lat0 = lat0_deg * degtorad
    lon0 = lon0_deg * degtorad
    clat,clon = np.cos(lats), np.cos(lons)
    slat,slon = np.sin(lats), np.sin(lons)
    delX = np.cos(lat0)*np.cos(lon0) - clat*clon
    delY = np.cos(lat0)*np.sin(lon0) - clat*slon
    delZ = np.sin(lat0) - slat;
    dist_sq = np.square(delX) + np.square(delY) + np.square(delZ)
    minindex_1d = dist_sq.argmin()  # 1D index of minimum element
    iy_min,ix_min = np.unravel_index(minindex_1d, lats.shape)
    return (iy_min,ix_min)

class File_Geom (object):
    def __init__(self, filename, groupname=None):
        nc = NetCDFFile(filename, 'r')
        if not groupname is None:
            grp = nc.groups[groupname]
        else:
            grp = nc
        self.nc = nc
        self.lon = grp.variables['longitude'][:,:]
        self.lat = grp.variables['latitude'][:,:]
        self.lon_bnds = grp.variables['longitude_bounds'][:,:,:]
        self.lat_bnds = grp.variables['latitude_bounds'][:,:,:]

    def __del__(self):
        self.nc.close()

    def get_patches (self, lon0, lat0, size):
        iy,ix = find_point (self.lon, self.lat, lon0, lat0)
        shape = self.lat.shape
        hw = int(size/2)
        patches = []
        idlist = []
        k = 0
        for i in range(max(0,iy-hw),min(iy+hw,shape[0]-1)):
            for j in range(max(0,ix-hw),min(ix+hw,shape[1]-1)):
                p_coords = np.stack((self.lon_bnds[i,j,:], self.lat_bnds[i,j,:]), axis=1)
                poly = Polygon(p_coords, True)
                patches.append (poly)
                idlist.append(k)
                k += 1
        return patches, idlist

def main():
    parser = argparse.ArgumentParser(description='plot polygon vertices')
    parser.add_argument('--format', default='pdf',
                        help="plot format [e.g. pdf, svgz]")
    parser.add_argument('--center', metavar=('LON','LAT',), default=None, nargs=2, type=float,
                        help="center lon,lat [deg]")
    parser.add_argument('--size', default=20, type=float,
                        help="plot width [pixels]")
    parser.add_argument('--group', default=None,
                        help="group containing pixel corners")
    parser.add_argument('--outfile', default=None,
                        help="output plot file")
    parser.add_argument('filename', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    lon0 = args.center[0]
    lat0 = args.center[1]

    if args.outfile is None:
        figfile = "pixels." + args.format
    else:
        figfile = args.outfile

    fig, ax = plt.subplots()
    ax.ticklabel_format (style='sci', scilimits=(-3,4), useMathText=True)

    ax.set_xlabel ('longitude [deg]')
    ax.set_ylabel ('latitude [deg]')

    for fn in args.filename:
        pixels = File_Geom (fn, args.group)
        patches, idlist = pixels.get_patches (lon0, lat0, args.size)
        p = PatchCollection (patches, edgecolor='k', alpha=0.4, linewidth=0.125)
        p.set_array (np.array(idlist))
        ax.add_collection(p)

    ax.autoscale_view()

    print("Saving file: {}".format(figfile))
    fig.savefig (figfile, bbox_inches='tight')

if __name__ == "__main__":
    main()
