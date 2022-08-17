#! /home/temposdpc/miniconda3/bin/python

import os
# Docs suggest this might improve performance when threading doesn't matter.
# I'm not seeing it.
os.environ["PYPROJ_GLOBAL_CONTEXT"]="ON"
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.backends.backend_pdf
import cartopy.crs as ccrs
import cartopy.img_transform as ctr
import cartopy.feature as cfeature
from netCDF4 import Dataset as NetCDFFile

#earth_mean_radius=6371008.8 # meters

class Var_Map (object):
    def __init__(self, var, keys):
        self.var = var
        self.keys = keys

class Sza_File (object):

    def __init__(self, filename, central_longitude, central_latitude):
        nc = NetCDFFile(filename, 'r')
        self.nc = nc
        self.lon = nc.variables['longitude'][:]
        self.lat = nc.variables['latitude'][:]
        # Plate Carree projection
        self.xpc = nc.variables['x'][:,:]
        self.ypc = nc.variables['y'][:,:]
        # day begin/end control points
        self.day_beg_point = nc.getncattr ('day_begin_ctrl_point')
        self.day_end_point = nc.getncattr ('day_end_ctrl_point')
        self.plan_id = nc.getncattr ('plan_id')

        # Map projections
        globe = None
        self.gdt = ccrs.Geodetic ()
        self.eqc = ccrs.PlateCarree (central_longitude=central_longitude, globe=globe)
        self.nsp = ccrs.NearsidePerspective(central_longitude=central_longitude, central_latitude=central_latitude,
                                           satellite_height=35785831,
                                           false_easting=0, false_northing=0)
    def __del__(self):
        self.nc.close()

    def read_var (self, name):
        var_ptr = self.nc.variables[name]
        fill = var_ptr.getncattr ('_FillValue')
        var = var_ptr[:]
        var[var == fill] = np.nan
        keys = {}
        keys["jd_utc_str"] = var_ptr.getncattr('julian_date_str')
        keys["solar_boresight_angle"] = var_ptr.getncattr('solar_boresight_angle')
        keys["scan_duration"] = var_ptr.getncattr('scan_duration')
        keys["num_repeats_cbm"] = var_ptr.getncattr('num_repeats_cbm')
        keys["box_lon"] = var_ptr.getncattr('box_lon')
        keys["box_lat"] = var_ptr.getncattr('box_lat')
        keys["maneuver_loss"] = var_ptr.getncattr('maneuver_loss')
        return Var_Map (var, keys)

    def plot_var (self, var):
        sza = np.ma.masked_invalid (var.var)
        xx = np.ma.masked_invalid (self.xpc)
        yy = np.ma.masked_invalid (self.ypc)
        sza.mask = np.logical_or (sza.mask, xx.mask, yy.mask)

        extent = (np.min(xx), np.max(xx), np.min(yy), np.max(yy))
        eqc = self.eqc

        plot_array = sza
        # Warp to a different projection for plotting
        plot_proj = self.nsp
        (plot_array, plot_extent) = ctr.warp_array (sza, plot_proj, source_proj=eqc)

        num_levels=9
        bounds = np.linspace(0, 90, num_levels+1)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

        cmap_name='plasma_r' #'hot_r'

        fig = plt.figure (figsize=(7,3.5), dpi=100)
        fig.subplots_adjust (left=0.1, right=0.85, top=0.95, bottom=0.05)
        ax = plt.axes(projection=plot_proj)
        # set axes extent and labels in native projection
        ax.set_extent(extent, crs=eqc)
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth=0.2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.2)
        ax.add_feature(cfeature.STATES, linewidth=0.2)
        g = ax.gridlines(draw_labels=True)
        g.top_labels = False
        g.right_labels = False
        g.xlabel_style["size"] = 6
        g.ylabel_style["size"] = 6

        transform_first=True
        filled_c = ax.contourf(xx, yy, sza, levels=bounds,
                   transform=eqc, transform_first=transform_first,
                   cmap=plt.get_cmap(cmap_name))

        ## Create colorbar axes (temporarily) anywhere
        cax = fig.add_axes([0,0,0.1,0.1])
        ## Find the location of the main plot axes
        posn = ax.get_position()
        ## Adjust the positioning of the colorbar,
        cax.set_position([posn.x0+posn.width+0.02, posn.y0, 0.02, posn.height])

        cb = plt.colorbar(filled_c, ticks=bounds, cax=cax)
        cb.ax.tick_params(labelsize=6)
        cb.ax.set_ylabel ('SZA [deg]', fontsize=8)

        ax.contour(xx, yy, sza, levels=filled_c.levels,
                   colors=['black'], linewidths=0.1,
                   transform=eqc, transform_first=transform_first)

        box_lon = np.ma.masked_invalid (var.keys["box_lon"])
        box_lat = np.ma.masked_invalid (var.keys["box_lat"])

        box, = ax.plot (box_lon, box_lat, transform=self.gdt, color='red', linewidth=1)
        ax.plot (self.day_beg_point[0], self.day_beg_point[1], marker='.', transform=self.gdt, color='red')
        pt, = ax.plot (self.day_end_point[0], self.day_end_point[1], marker='.', linestyle='None', transform=self.gdt, color='red')

        ax.legend ([box, pt], ['Scan region', 'Control point'], loc='lower left',
                   bbox_to_anchor=(0.79, -0.28), borderaxespad=0.0, fontsize=8)

        if var.keys["maneuver_loss"] > 0.0:
            maneuver_truncation = "M:%0.1f sec" % (var.keys["maneuver_loss"])
        else:
            maneuver_truncation = ""

        ax.set_title ('%s  SBA:%0.1f deg %s' % (var.keys["jd_utc_str"],
                      var.keys["solar_boresight_angle"], maneuver_truncation), fontsize=8, loc='left')

        if var.keys["num_repeats_cbm"] == 0:
            ax.set_title ('%0.1f sec' % (var.keys["scan_duration"]), fontsize=8, loc='right')
        else:
            ax.set_title ('CBM:%d, %0.1f sec' % (var.keys["num_repeats_cbm"], var.keys["scan_duration"]),
                          fontsize=8, loc='right')

        ax.text (0.0, -0.2, r'plan_id = {}'.format(self.plan_id), fontsize=8, transform=ax.transAxes)

        return fig

def find_vars_for_date (nc, date):
    var_names = [key for key in nc.variables.keys() if key.startswith("sza_")]
    select_vars = []
    for name in var_names:
        var = nc.variables[name]
        date_str = var.getncattr('julian_date_str')
        if date_str.find(date) == 0:
            select_vars.append(name)
    return select_vars

def select_maneuver_affected (nc, var_names):
    select_vars = []
    for name in var_names:
        var = nc.variables[name]
        if var.getncattr('maneuver_loss') > 0.0:
            select_vars.append(name)
    return select_vars

def main():
    parser = argparse.ArgumentParser(description='Generate SZA plots.')
    parser.add_argument ('--date', metavar='YYYY-MM-DD', default=None,
                         help="Select plots by date")
    parser.add_argument ('--maneuver', action='store_true',
                         help="Plot only scans affected by maneuvers", )
    parser.add_argument ('--select', metavar='N', default=None, nargs="*", type=int,
                         help="Select plots by number")
    parser.add_argument ('--output', metavar='FILE', default="sza.pdf",
                         help="Output plot file name")
    parser.add_argument ('--central-latitude', metavar='LON', default=5.0,
                         help="Plotted nearside projection central latitude [deg]")
    parser.add_argument ('--central-longitude', metavar='LON', default=-90.0,
                         help="Input SZA map Plate Carree projection central longitude [deg]")
    parser.add_argument ('szafile', help="Path to SZA file (plan output)")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    s = Sza_File (args.szafile, args.central_longitude, args.central_latitude)

    if args.date is not None:
        var_name_list = find_vars_for_date (s.nc, args.date)
    elif args.select is not None:
        var_name_list = ['sza_%02d' % (num) for num in args.select]
    else:
        var_name_list = [key for key in s.nc.variables.keys() if key.startswith("sza_")]

    if args.maneuver:
        var_name_list = select_maneuver_affected (s.nc, var_name_list)

    num_vars = len(var_name_list)

    if num_vars == 0:
        print('No plots selected')
        sys.exit(0)

    print('Writing {} plots to: {}'.format(num_vars, args.output))
    pdf = matplotlib.backends.backend_pdf.PdfPages(args.output)

    k=1
    for var_name in var_name_list:
        print('Plotting {}  [{}/{}]'.format(var_name, k, num_vars))
        var = s.read_var (var_name)
        fig = s.plot_var(var)
        pdf.savefig (fig, bbox_inches='tight')
        plt.close(fig)
        k += 1
    pdf.close()

if __name__ == "__main__":
    main()

