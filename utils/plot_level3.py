#! /usr/bin/env python3

import sys
import glob
import os
import argparse

import matplotlib
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
    def __init__(self, lon, lat, var, units):
        self.lon = lon
        self.lat = lat
        self.var = var
        self.units = units

class Var_Map_Config (object):
    def __init__(self, varpath, min, max):
        self.name = os.path.basename(varpath)
        self.min = min
        self.max = max
        self.vargroup = os.path.dirname(varpath)

# Note that to plot the granule images without overlap,
# pcolormesh needs the xi,yi coordinate array dimensions
# to be one larger than the image array dimensions so that
# the outer pixel edge locations are known.
def read_var (filename, var_config, layer):
    nc = NetCDFFile(filename)
    lon = nc['/longitude'][:]
    lat = nc['/latitude'][:]
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
    return Var_Map (lon, lat, var, units)

def main():
    parser = argparse.ArgumentParser(description='plot science data')
    parser.add_argument('--outfile', help="path to output file", default=None)
    parser.add_argument('--varpath', help="path to variable in file", default=None)
    parser.add_argument('--varmin', help="min plot value", type=float)
    parser.add_argument('--varmax', help="max plot value", type=float)
    parser.add_argument('--layer', help="", default=None)
    parser.add_argument('--label', help="", default=None)
    parser.add_argument('filename', nargs=1)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.label == None:
        extra_label = ""
    else:
        extra_label = args.label

    plot_var = (args.varpath is not None)

    filename = args.filename[0]

    outfile = args.outfile
    if outfile is not None:
        matplotlib.use('Agg')

    # parse OGS boundary polygon attribute
    with NetCDFFile(filename) as nc:
        geospatial_bounds = nc.getncattr("geospatial_bounds")
    bdry_pts = geospatial_bounds.lstrip('POLYGON((').rstrip('))').split(',')
    bdry_lon = [float(p.split(' ')[1]) for p in bdry_pts]
    bdry_lat = [float(p.split(' ')[0]) for p in bdry_pts]

    if plot_var:
        var_config = Var_Map_Config (args.varpath, args.varmin, args.varmax)

    fig = plt.figure()
    ax = plt.subplot (1,1,1, projection=ccrs.Miller())
    if plot_var:
        ax.set_extent ([-155, -30, 15, 62], ccrs.PlateCarree())
    else:
        ax.set_extent ([min(bdry_lon), max(bdry_lon), min(bdry_lat), max(bdry_lat)], ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.5)
    gl.top_labels=False
    gl.right_labels=False
    gl.xlabel_style = {'size':6}
    gl.ylabel_style = {'size':6}
    ax.add_feature (cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature (cartopy.feature.BORDERS, linewidth=0.5)
    ax.add_feature (cartopy.feature.STATES, linewidth=0.5)

    cmap = plt.get_cmap('jet')

    if plot_var:
        vm = read_var (filename, var_config, args.layer)
        lons, lats = np.meshgrid (vm.lon, vm.lat)
        cs = ax.pcolormesh (lons, lats, vm.var, cmap=cmap, rasterized=True,
                            vmin=var_config.min, vmax=var_config.max,
                            transform=ccrs.PlateCarree())

        sm = plt.cm.ScalarMappable (norm=cs.norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar (sm, ax=ax, orientation='horizontal', pad=0.05, aspect=50, format='%.2e')
        cbar.set_label(vm.units, size=6)
        cbar.ax.tick_params(labelsize=6)

    ax.plot (bdry_lon, bdry_lat, color='g', lw=0.5, transform=ccrs.PlateCarree())
    plt.title ("{}".format(os.path.basename(filename)), fontsize='x-small')
    if plot_var:
        plt.suptitle ("{}{}".format(var_config.name, extra_label), y=0.825)

    if outfile is None:
        plt.show()
    else:
        print('Creating plot file: {}'.format(outfile))
        plt.savefig(outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

