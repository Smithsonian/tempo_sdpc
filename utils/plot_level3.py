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
    parser.add_argument('--outfile', help="path to output file")
    parser.add_argument('--varpath', help="path to variable in file")
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

    filename = args.filename[0]

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

    plt.suptitle ("{}{}".format(var_config.name, extra_label), y=0.73)
    plt.title ("{}".format(os.path.basename(filename)), fontsize='x-small')

    #plt.show()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

