#! /usr/bin/env python3

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
    def __init__(self, var, jd_utc, jd_utc_str, solar_boresight_angle, scan_duration, num_repeats, num_repeats_cbm, start_pos, scan_angle, box_lon, box_lat):
        self.var = var
        self.jd_utc = jd_utc
        self.jd_utc_str = jd_utc_str
        self.solar_boresight_angle = solar_boresight_angle
        self.scan_duration = scan_duration
        self.num_repeats = num_repeats
        self.num_repeats_cbm = num_repeats_cbm
        self.start_pos = start_pos
        self.scan_angle = scan_angle
        self.box_lon = box_lon
        self.box_lat = box_lat

class Sza_File (object):

    def __init__(self, filename):
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

        # Map projections
        central_longitude=-90.0
        globe = None
        self.gdt = ccrs.Geodetic ()
        self.eqc = ccrs.PlateCarree (central_longitude=central_longitude, globe=globe)
        self.nsp = ccrs.NearsidePerspective(central_longitude=central_longitude, central_latitude=0.0,
                                           satellite_height=35785831,
                                           false_easting=0, false_northing=0)
    def __del__(self):
        self.nc.close()

    def read_var (self, name):
        var_ptr = self.nc.variables[name]
        fill = var_ptr.getncattr ('_FillValue')
        var = var_ptr[:]
        var[var == fill] = np.nan
        jd_utc = var_ptr.getncattr('julian_date')
        jd_utc_str = var_ptr.getncattr('julian_date_str')
        solar_boresight_angle = var_ptr.getncattr('solar_boresight_angle')
        scan_duration = var_ptr.getncattr('scan_duration')
        num_repeats = var_ptr.getncattr('num_repeats')
        num_repeats_cbm = var_ptr.getncattr('num_repeats_cbm')
        start_pos = var_ptr.getncattr('start_pos')
        scan_angle = var_ptr.getncattr('scan_angle_rad')
        box_lon = var_ptr.getncattr('box_lon')
        box_lat = var_ptr.getncattr('box_lat')
        return Var_Map (var, jd_utc, jd_utc_str, solar_boresight_angle, scan_duration, num_repeats, num_repeats_cbm, start_pos, scan_angle, box_lon, box_lat)

    def plot_var (self, var):
        sza = var.var
        xx = self.xpc
        yy = self.ypc
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

        box_lon = np.ma.masked_invalid (var.box_lon)
        box_lat = np.ma.masked_invalid (var.box_lat)

        ax.plot (box_lon, box_lat, transform=self.gdt, color='red', linewidth=1)
        ax.plot (self.day_beg_point[0], self.day_beg_point[1], marker='.', transform=self.gdt, color='red')
        ax.plot (self.day_end_point[0], self.day_end_point[1], marker='.', transform=self.gdt, color='red')

        ax.set_title ('%s  sba:%0.1f deg' % (var.jd_utc_str, var.solar_boresight_angle), fontsize=8, loc='left')
        if var.num_repeats_cbm == 0:
            ax.set_title ('{} sec'.format(var.scan_duration), fontsize=8, loc='right')
        else:
            ax.set_title ('cbm:{}, {} sec'.format(var.num_repeats_cbm, var.scan_duration),
                          fontsize=8, loc='right')
        return fig

def main():
    parser = argparse.ArgumentParser(description='Generate SZA plots.')
    parser.add_argument ('--output', metavar='FILE', default="sza.pdf",
                         help="Output plot file name")
    parser.add_argument('--select', metavar='LIST', default=None, nargs="*", type=int,
                        help="Selected plot numbers")
    parser.add_argument ('szafile', help="Path to SZA file (plan output)")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    s = Sza_File (args.szafile)

    pdf = matplotlib.backends.backend_pdf.PdfPages(args.output)

    if args.select is not None:
        var_name_list = ['sza_%02d' % (num) for num in args.select]
    else:
        var_name_list = [key for key in s.nc.variables.keys() if key.startswith("sza_")]

    for var_name in var_name_list:
        print('Plotting {}'.format(var_name))
        var = s.read_var (var_name)
        fig = s.plot_var(var)
        pdf.savefig (fig)
        plt.close(fig)
    pdf.close()

if __name__ == "__main__":
    main()

