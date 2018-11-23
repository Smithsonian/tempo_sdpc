#! /usr/bin/python3

import sys, os
import argparse
parser = argparse.ArgumentParser(description='Make a plot.')
parser.add_argument('filename', help="netCDF data file name")
args = parser.parse_args()

nc_filename = args.filename
plt_filename = os.path.splitext(nc_filename)[0] + ".pdf"

import numpy as np
import numpy.ma as ma
import matplotlib
# matplotlib.use('Agg')
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile

class Grid_Map (object):
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat

class Var_Map (object):
    def __init__(self, var, jd_utc, jd_utc_str, scan_duration, num_repeats, start_pos, scan_angle):
        self.var = var
        self.jd_utc = jd_utc
        self.jd_utc_str = jd_utc_str
        self.scan_duration = scan_duration
        self.num_repeats = num_repeats
        self.start_pos = start_pos
        self.scan_angle = scan_angle

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
    scan_duration = var_ptr.getncattr('scan_duration')
    num_repeats = var_ptr.getncattr('num_repeats')
    start_pos = var_ptr.getncattr('start_pos')
    scan_angle = var_ptr.getncattr('scan_angle_rad')
    return Var_Map (var, jd_utc, jd_utc_str, scan_duration, num_repeats, start_pos, scan_angle)

def plot_var_map (m, xi, yi, var, var_config):
    clevs = np.arange(var_config.min, var_config.max, 10)
    m.contour (xi,yi,var.var, [90.0], linewidths=2)
    # For greyscale, try cmap='bone_r'
    # For color, try cmap = 'bwr_r', 'YlOrBr', or 'hsv'
    cs = m.contourf(xi,yi,var.var, clevs,
                    cmap='plasma_r', alpha=0.625, extend='both')

    # scan start line
    (x0, y0) = m(var.start_pos[0], var.start_pos[1]) # (lon,lat) -> (x,y)
    ones_i = np.ones(len(yi))
    scan_reg_color='white'
    scan_reg_linewidth=2.5
    plt.plot (x0 * ones_i, yi, color=scan_reg_color, linewidth=scan_reg_linewidth)
    (xc, yc) = m(-100.0, 36.0)
    plt.plot (x0, yc, marker=8, color=scan_reg_color, markersize=15, markeredgewidth=scan_reg_linewidth, fillstyle='none')

    # scan end line
    geo_altitude = 35785831.0  # meters
    x1 = x0 + var.scan_angle * geo_altitude;
    plt.plot (x1 * ones_i, yi, color=scan_reg_color, linewidth=scan_reg_linewidth, linestyle='--')

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
    #x_lab = [-124, -126, -130]
    x_lab = [-128.75, -132, -137.5, -148.5]
    for i in range(len(x_lab)):
        plt.annotate(r"\it %g N" % (parallels[i]),
                     xy=m(x_lab[i],parallels[i]),xycoords='data',
                     horizontalalignment='right',verticalalignment='center',
                     fontsize=15,annotation_clip=False)
    #######
    #parallels = np.arange(0.0,90,10.)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    #m.drawmeridians(meridians) # cannot label meridians on 'ortho' projections
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15, fmt=(lambda x: (r"\it %d W" % (abs(x-360)))))

def plot_var (nc, m, xi, yi, var_config):
    config_map (m)
    var = read_var (nc, var_config)
    cs = plot_var_map (m, xi, yi, var, var_config)
    cbar = m.colorbar(cs,location='bottom',pad="10%")
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(var_config.cbar_label, size=15)
    plt.title(var.jd_utc_str, loc='left', size=15)
    plt.title(r"%d$\times$ %0.1f min" % (var.num_repeats, var.scan_duration/60.0), loc='right', size=15)

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

nc = NetCDFFile(nc_filename)
grid = read_grid (nc)
xi, yi = m(grid.lon, grid.lat)

var_names = list(nc.variables.keys())
sza_var_indices = [i for i,item in enumerate(var_names) if "sza_" in item]

pdf = matplotlib.backends.backend_pdf.PdfPages(plt_filename)

for k in sza_var_indices:
    fig = plt.figure(1)
    fig.set_size_inches (9.0, 6.5)
    fig.set_dpi (100)
    fig.suptitle ('Scan endpoints', fontsize=20, x=0.515)
    print(var_names[k])
    config = Var_Map_Config (var_names[k], 'SZA [deg]', 50, 140, 'xxx')
    plot_var (nc, m, xi, yi, config)
    pdf.savefig(fig)
    plt.close(fig)

pdf.close()
nc.close()

