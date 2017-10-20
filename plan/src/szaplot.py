#! /usr/bin/python3

import sys, os
import argparse
parser = argparse.ArgumentParser(description='Make a plot.')
parser.add_argument('filename', help="netCDF data file name")
parser.add_argument('-x', '--ext', dest='ext', default="pdf",
                      help="plot type, EXT=pdf | png | ps | ...")
args = parser.parse_args()

nc_filename = args.filename
plt_filename = os.path.splitext(nc_filename)[0] + "." + args.ext

import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile

class Grid_Map (object):
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat

class Var_Map (object):
    def __init__(self, var, jd_utc, jd_utc_str):
        self.var = var
        self.jd_utc = jd_utc
        self.jd_utc_str = jd_utc_str

class Var_Map_Config (object):
    def __init__(self, name, cbar_label, min, max, title):
        self.name = name
        self.cbar_label = cbar_label
        self.title = title
        self.min = min
        self.max = max

def read_grid (nc):
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    return Grid_Map (lon, lat)

def read_var (nc, var_config):
    var_ptr = nc.variables[var_config.name]
    var = var_ptr[:]
    jd_utc = var_ptr.getncattr('julian_date')
    jd_utc_str = var_ptr.getncattr('julian_date_str')
    return Var_Map (var, jd_utc, jd_utc_str)

def plot_var_map (m, grid, var, var_config):
    xi, yi = m(grid.lon, grid.lat)
    clevs = np.arange(var_config.min, var_config.max, 10)
    # For greyscale, try cmap='bone_r'
    # For color, try cmap = 'bwr_r', 'YlOrBr', or 'hsv'
    cs = m.contourf(xi,yi,var.var, clevs, cmap='bwr_r', extend='both')
    return cs

def config_map (m):
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines(linewidth=0.25)
    m.drawstates(linewidth=0.25)
    m.drawcountries(linewidth=0.25)
    # draw parallels.
    parallels = np.arange(20.0,60.0,10.)
    m.drawparallels(parallels)
    #### for 'geos', 'ortho' projections, label parallels manually (UGLY!)
    ## x_lab = [-120, -120]
    ## for i in range(len(x_lab)):
    ##     plt.annotate(np.str(parallels[i]),xy=m(x_lab[i],parallels[i]),xycoords='data',
    ##                  horizontalalignment='right')
    #######
    #parallels = np.arange(0.0,90,10.)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    #m.drawmeridians(meridians) # cannot label meridians on 'ortho' projections
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

def plot_var (nc, m, grid, var_config):
    config_map (m)
    var = read_var (nc, var_config)
    cs = plot_var_map (m, grid, var, var_config)
    cbar = m.colorbar(cs,location='bottom',pad="10%")
    cbar.set_label(var_config.cbar_label)
    plt.title('{}\n{}'.format (var_config.title, var.jd_utc_str))

# Create Basemap instance.
# 'geos' seems to be the best projection for this application, but here
# are some other things I tried:
#
# m = Basemap(projection='stere',lon_0=-95.0,lat_0=36.0,lat_ts=42.0,resolution='l',
#              width=7500000,height=7500000)
# m = Basemap(projection='ortho',lon_0=-100,lat_0=36.0,resolution='l')
# m = Basemap(projection='nsper',satellite_height=4000.0e3,
#             lon_0=-95.0,lat_0=36.0,resolution='l')
# m = Basemap(projection='aea',lon_0=-100.0,lat_0=42.0,lat_ts=42.0,resolution='l',
#             llcrnrlon=-125,llcrnrlat=8,urcrnrlon=-30,urcrnrlat=60)

# For the 'geos' projection, use m1 to define the boundaries:
lon_0=-100.0
m1 = Basemap(projection='geos',lon_0=lon_0,resolution=None)
px = m1.urcrnrx * 0.3
py = m1.urcrnry * 0.48
mx = m1.urcrnrx * (-0.25)
my = m1.urcrnry * 0.15
m  = Basemap(projection='geos',lon_0=lon_0,resolution='l',\
     llcrnrx=mx,llcrnry=my,urcrnrx=px,urcrnry=py)

beg_config = Var_Map_Config ('sza_beg', 'SZA [deg]', 50, 140, 'First scan, start')
end_config = Var_Map_Config ('sza_end', 'SZA [deg]', 50, 140, 'Last scan, end')

nc = NetCDFFile(nc_filename)
grid = read_grid (nc)

fig = plt.figure(1)
fig.set_size_inches (9.0, 6.5)
fig.set_dpi (200)
fig.suptitle ('Solar zenith angle', fontsize=20, y=0.78)
ax = plt.subplot(121)
plot_var (nc, m, grid, end_config)
ax = plt.subplot(122)
plot_var (nc, m, grid, beg_config)

nc.close()

#plt.show()
fig.savefig (plt_filename, orientation='landscape',
             bbox_inches='tight', pad_inches=0.1)
plt.close(fig)

