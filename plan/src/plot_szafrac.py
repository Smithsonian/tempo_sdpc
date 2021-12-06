#! /usr/bin/python3

import sys, os
import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset as NetCDFFile

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

class Var_Map (object):
    def __init__(self, var, jd_utc, jd_utc_str, scan_duration, num_repeats, start_pos, scan_angle):
        self.var = var
        self.jd_utc = jd_utc
        self.jd_utc_str = jd_utc_str
        self.scan_duration = scan_duration
        self.num_repeats = num_repeats
        self.start_pos = start_pos
        self.scan_angle = scan_angle

def read_var (nc, name):
    var_ptr = nc.variables[name]
    var = var_ptr[:]
    jd_utc = var_ptr.getncattr('julian_date')
    jd_utc_str = var_ptr.getncattr('julian_date_str')
    scan_duration = var_ptr.getncattr('scan_duration')
    num_repeats = var_ptr.getncattr('num_repeats')
    start_pos = var_ptr.getncattr('start_pos')
    scan_angle = var_ptr.getncattr('scan_angle_rad')
    return Var_Map (var, jd_utc, jd_utc_str, scan_duration, num_repeats, start_pos, scan_angle)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Plot SZA distribution at the start of each scan')
    parser.add_argument('infile', help="netCDF data file name")
    parser.add_argument('--outfile', help="plot file name")
    if len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    nc_filename = args.infile

    nc = NetCDFFile(nc_filename, 'r')
    var_names = list(nc.variables.keys())
    sza_vars = [var_names[i] for i,item in enumerate(var_names) if "sza_" in item]

    pdf = PdfPages (args.outfile)
    fig, ax1 = plt.subplots (1,1)
    fig = plt.figure(1)

    ax1.set_ylim (0, 1.0)
    ax1.set_xlabel ('mirror x [mrad]')
    ax1.set_ylabel ('f$_y$ (SZA\,$<70$\,deg)')
    ax1.set_title ('SZA @ scan start')

    # convert radian to milliradian
    x = nc.variables['x'][:] * 1.e3
    # factor of 2 between mirror rotation angle coordinate
    # and angular step in FOR for a reflected ray
    x /= 2
    for var_name in sza_vars:
        var = read_var (nc, var_name)
        lit = np.where(var.var < 70.0, 1, 0)
        y = np.sum(lit, axis=0)
        ax1.plot(x, y / var.var.shape[0]);

    pdf.savefig(fig)
    pdf.close()

    nc.close()

if __name__ == '__main__':
    main()

