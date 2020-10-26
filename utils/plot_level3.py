#! /usr/bin/env python

import sys
import glob
import os
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile
import numpy as np

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

def config_map (m):
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    # draw parallels.
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

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

    # create Basemap instance.
    # m = Basemap(projection='stere',lon_0=-95.0,lat_0=48.0,lat_ts=42.0,resolution='l',
    #             width=7500000,height=7500000)
    # m = Basemap(projection='ortho',lon_0=-95,lat_0=42,resolution='l')
    # m = Basemap(projection='aea',lon_0=-90.0,lat_0=42.0,lat_ts=42.0,resolution='l',
    #             llcrnrlon=-125,llcrnrlat=8,urcrnrlon=-30,urcrnrlat=60)

    m = Basemap(projection='mill',lon_0=-90.0,lat_0=40.0,resolution='l',
    llcrnrlon=-155,llcrnrlat=15,urcrnrlon=-30,urcrnrlat=65)

    config_map (m)

    var_config = Var_Map_Config (args.varpath, args.varmin, args.varmax)

    cmap = plt.get_cmap('jet')

    vm = read_var (filename, var_config, args.layer)
    lons, lats = np.meshgrid (vm.lon, vm.lat)
    xi, yi = m(lons, lats)
    cs = m.pcolormesh (xi, yi, vm.var, cmap=cmap, rasterized=True,
                       vmin=var_config.min, vmax=var_config.max)

    cbar = m.colorbar(cs,location='right',pad="4%")
    cbar.set_label(vm.units)
    plt.suptitle ("{}{}".format(var_config.name, extra_label), y=0.8)
    plt.title ("{}".format(os.path.basename(filename)), fontsize='x-small')

    #plt.show()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

