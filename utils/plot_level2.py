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
    def __init__(self, lon_bnds, lat_bnds, var, units):
        self.lon_bnds = lon_bnds
        self.lat_bnds = lat_bnds
        self.var = var
        self.units = units

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
    return Var_Map (lon_bnds, lat_bnds, var, units)

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
def plot_var_map (m, vm, var_config, cmap):
    for ix in range(vm.lon_bnds.shape[0]):
        lons = select_scan_step (vm.lon_bnds, ix)
        lats = select_scan_step (vm.lat_bnds, ix)
        xi, yi = m(lons, lats)
        var = np.transpose(vm.var[[ix],:])
        # If not rasterized, the output plot would be enormous.
        cs = m.pcolormesh(xi, yi, var, cmap=cmap, rasterized=True,
                          vmin=var_config.min, vmax=var_config.max)
    return cs

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
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.label == None:
        extra_label = ""
    else:
        extra_label = args.label

    # create Basemap instance.
    # m = Basemap(projection='stere',lon_0=-95.0,lat_0=48.0,lat_ts=42.0,resolution='l',
    #             width=7500000,height=7500000)
    # m = Basemap(projection='ortho',lon_0=-95,lat_0=42,resolution='l')
    # m = Basemap(projection='aea',lon_0=-90.0,lat_0=42.0,lat_ts=42.0,resolution='l',
    #             llcrnrlon=-125,llcrnrlat=8,urcrnrlon=-30,urcrnrlat=60)

    m = Basemap(projection='mill',lon_0=-90.0,lat_0=40.0,resolution='l',
    llcrnrlon=-155,llcrnrlat=15,urcrnrlon=-30,urcrnrlat=65)

    config_map (m)

    filenames = args.filenames
    filenames.sort()

    var_config = Var_Map_Config (args.varpath, args.varmin, args.varmax)

    cmap = plt.get_cmap('jet')

    for f in filenames:
        print('reading {}'.format(f))
        vm = read_var (f, var_config, args.layer)
        cs = plot_var_map (m, vm, var_config, cmap)

    cbar = m.colorbar(cs,location='right',pad="5%")
    cbar.set_label(vm.units)
    plt.title ("{}, ...\n{}".format(os.path.basename(filenames[0]),
                                    os.path.basename(filenames[-1])),
                                    fontsize='x-small')
    plt.suptitle("{}{}".format(var_config.name, extra_label), y=0.85)

    #plt.show()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == '__main__':
    main()

