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

class Var_Map (object):
    def __init__(self, lon_bnds, lat_bnds, var, units, bdry_lon, bdry_lat):
        self.lon_bnds = lon_bnds
        self.lat_bnds = lat_bnds
        self.var = var
        self.units = units
        self.bdry_lon = bdry_lon
        self.bdry_lat = bdry_lat

class Var_Map_Config (object):
    def __init__(self, varpath, min, max, geogrp='geolocation'):
        self.name = os.path.basename(varpath)
        self.min = min
        self.max = max
        self.vargroup = os.path.dirname(varpath)
        self.geogroup = geogrp

# Note that to plot the granule images without overlap,
# pcolormesh needs the xi,yi coordinate array dimensions
# to be one larger than the image array dimensions so that
# the outer pixel edge locations are known.
def read_var (filename, var_config, layer):
    nc = NetCDFFile(filename)
    geospatial_bounds = nc.getncattr("geospatial_bounds")
    bdry_pts = geospatial_bounds.lstrip('POLYGON((').rstrip('))').split(',')
    bdry_lon = [float(p.split(' ')[1]) for p in bdry_pts]
    bdry_lat = [float(p.split(' ')[0]) for p in bdry_pts]
    geogrp = nc.groups[var_config.geogroup];
    lon_bnds = geogrp.variables['longitude_bounds'][:]  # NE, NW, SW, SE
    lat_bnds = geogrp.variables['latitude_bounds'][:]
    vargrp = nc[var_config.vargroup]
    var_obj = vargrp.variables[var_config.name]
    var = var_obj[:]
    if layer:
        var = var[:,:,layer]
    if 'units' in var_obj.__dict__:
        units = var_obj.units
    else:
        units = ''
    nc.close()
    if var_config.min is None:
        var_config.min = var.min()
    if var_config.max is None:
        var_config.max = var.max()
    return Var_Map (lon_bnds, lat_bnds, var, units, bdry_lon, bdry_lat)

def select_scan_step (packed_corners, ix):
    # pixel corner packing sequence,
    # e.g. [NE, NW, SW, SE] => ne=0, nw=1, sw=2, se=3
    ne=0; nw=1; sw=2; se=3;
    # NW corner of the image is 0,0
    # SE corner of the image is -1,-1
    sides = np.dstack ((packed_corners[ix ,:,nw], packed_corners[ix, :,ne])).squeeze()
    bottom = packed_corners[ix,-1,[sw,se]]
    result = np.vstack ((sides, bottom))
    return result

# It would be nice if pcolormesh could mask points with invalid coordinates.
# Apparently, that's not supported.
def plot_var_map (ax, vm, var_config, cmap):
    for ix in range(vm.lon_bnds.shape[0]):
        lons = select_scan_step (vm.lon_bnds, ix)
        lats = select_scan_step (vm.lat_bnds, ix)
        var = np.transpose(vm.var[[ix],:])
        # If not rasterized, the output plot would be enormous.
        cs = ax.pcolormesh(lons, lats, var, cmap=cmap, rasterized=True,
                           vmin=var_config.min, vmax=var_config.max,
                           transform=ccrs.PlateCarree())
        ax.plot (vm.bdry_lon, vm.bdry_lat, color='g', lw=0.5, transform=ccrs.PlateCarree())
    return cs

def main():
    parser = argparse.ArgumentParser(description='plot science data')
    parser.add_argument('--outfile', help="path to output file")
    parser.add_argument('--varpath', help="path to variable in file")
    parser.add_argument('--varmin', help="min plot value", type=float)
    parser.add_argument('--varmax', help="max plot value", type=float)
    parser.add_argument('--layer', help="", default=None, type=int)
    parser.add_argument('--label', help="", default=None)
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.label == None:
        extra_label = ""
    else:
        extra_label = args.label

    filenames = args.filenames
    filenames.sort()

    var_config = Var_Map_Config (args.varpath, args.varmin, args.varmax)

    fig = plt.figure()
    ax = plt.subplot (1,1,1, projection=ccrs.Miller())
    ax.set_extent ([-155, -30, 15, 55], ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5)
    gl.top_labels=False
    gl.right_labels=False
    gl.xlabel_style = {'size':6}
    gl.ylabel_style = {'size':6}
    ax.add_feature (cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature (cartopy.feature.BORDERS, linewidth=0.5)
    ax.add_feature (cartopy.feature.STATES, linewidth=0.5)

    cmap = plt.get_cmap('jet')

    for f in filenames:
        print('reading {}'.format(f))
        vm = read_var (f, var_config, args.layer)
        cs = plot_var_map (ax, vm, var_config, cmap)

    sm = plt.cm.ScalarMappable (norm=cs.norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar (sm, ax=ax, orientation='horizontal', pad=0.05, aspect=50, format='%.2e')
    cbar.set_label(vm.units, size=6)
    cbar.ax.tick_params(labelsize=6)

    plt.title ("{}, ...\n{}".format(os.path.basename(filenames[0]),
                                    os.path.basename(filenames[-1])),
                                    fontsize='x-small')
    plt.suptitle("{}{}".format(var_config.name, extra_label), y=0.75)

    #plt.show()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

