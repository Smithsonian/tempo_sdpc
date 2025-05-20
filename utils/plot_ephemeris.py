#! /usr/bin/env python3

import sys
import csv
import numpy as np
import re
import time
import argparse

from netCDF4 import Dataset

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import gridspec

import matplotlib.ticker as ticker
import matplotlib.dates as dates

plt.rc('text', usetex=True)

def append_fields (a, b):
    if len(a.keys()) == 0:
        return b
    for key, value in a.items():
        a[key] = np.concatenate ((value, b[key]))
    return a

def sort_unique_ephem (eph):
    t = np.asarray (eph["t"])
    _t, i = np.unique (t, return_index=True)
    for key, value in eph.items():
        v = np.asarray(value)
        eph[key] = v[i]
    return eph

def read_dop_ephem (path):
    eph = {}
    with Dataset (path, "r") as nc:
        grp = nc.groups['ephemeris']
        eph["t"] = grp['dop_time'][:]
        eph["x"] = grp['dop_x'][:]
        eph["y"] = grp['dop_y'][:]
        eph["z"] = grp['dop_z'][:]
        eph["vx"] = grp['dop_vx'][:]
        eph["vy"] = grp['dop_vy'][:]
        eph["vz"] = grp['dop_vz'][:]
    return eph

def read_gpsr_ephem (path):
    eph = {}
    with Dataset (path, "r") as nc:
        grp = nc.groups['ephemeris']
        eph["t"] = grp['gpsr_time'][:]
        eph["x"] = grp['gpsr_x'][:]
        eph["y"] = grp['gpsr_y'][:]
        eph["z"] = grp['gpsr_z'][:]
        eph["vx"] = grp['gpsr_vx'][:]
        eph["vy"] = grp['gpsr_vy'][:]
        eph["vz"] = grp['gpsr_vz'][:]
        eph["mth"] = grp['gpsr_navsoln_mth'][:]
    return eph

def time_select_ephem (eph, start, end):
    if start is None:
        start = -np.Inf
    if end is None:
        end = np.Inf
    t = eph["t"]
    k = np.nonzero (np.logical_and(start <= t, t <= end))[0]
    for key, value in eph.items():
        eph[key] = value[k]

def main():
    parser = argparse.ArgumentParser(description='plot ECEF ephemeris comparison')
    parser.add_argument('--outfile', help="path to output PDF file")
    parser.add_argument('--start', type=float, default=None,
                        help="start time [sec since TEMPO epoch]")
    parser.add_argument('--end', type=float, default=None,
                        help="end time [sec since TEMPO epoch]")
    parser.add_argument('--dx', type=float, default=20.0, help="+/-X plot limits [km]")
    parser.add_argument('--dy', type=float, default=5.0, help="+/-Y plot limits [km]")
    parser.add_argument('--dz', type=float, default=35.0, help="+/-Z plot limits [km]")
    parser.add_argument('--dd', type=float, default=10.0, help="+/-Z plot limits [km]")
    parser.add_argument('files', nargs=argparse.REMAINDER,
                        help="path to netcdf4/HDF5 ephemeris data file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    outfile = args.outfile

    dop = {}
    gps = {}
    for eph_file in args.files:
        dop = append_fields (dop, read_dop_ephem (eph_file))
        gps = append_fields (gps, read_gpsr_ephem (eph_file))

    dop = sort_unique_ephem (dop)
    gps = sort_unique_ephem (gps)

    time_select_ephem (dop, args.start, args.end)
    time_select_ephem (gps, args.start, args.end)

    # hack to convert tempo timestamp to timet
    delta_taix_timet = 315964782.0;
    tt_dop = dop["t"] + delta_taix_timet
    tt_gps = gps["t"] + delta_taix_timet

    lw_std = 0.25

    with PdfPages (outfile) as pp:
        fig = plt.figure(figsize=(9, 6.5))

        gs = gridspec.GridSpec(4, 1, wspace=0.0, hspace=0.0, top=0.95, bottom=0.15, left=0.17, right=0.845)
        xax = plt.subplot(gs[0])
        yax = plt.subplot(gs[1], sharex=xax)
        zax = plt.subplot(gs[2], sharex=xax)
        dlt = plt.subplot(gs[3], sharex=xax)

        xax.tick_params(axis="x", direction='in', length=4, labelbottom=False)
        yax.tick_params(axis="x", direction='in', length=4, labelbottom=False)
        yax.ticklabel_format(axis="y", style='plain', useOffset=False)
        zax.tick_params(axis="x", direction='in', length=4, labelbottom=False)
        dlt.tick_params(axis="x", direction='in', length=4)

        x_med = np.median(gps["x"])
        y_med = np.median(gps["y"])
        z_med = np.median(gps["z"])
        xax.set_ylim (x_med - args.dx, x_med + args.dx)
        yax.set_ylim (y_med - args.dy, y_med + args.dy)
        zax.set_ylim (z_med - args.dz, z_med + args.dz)

        plt.xlim (tt_gps[0], tt_gps[-1])
        xax.plot (tt_gps, gps["x"], color='k', lw=lw_std, label='GPSR')
        yax.plot (tt_gps, gps["y"], color='k', lw=lw_std)
        zax.plot (tt_gps, gps["z"], color='k', lw=lw_std)

        xax.plot (tt_dop, dop["x"], color='r', lw=lw_std, label='DOP')
        yax.plot (tt_dop, dop["y"], color='r', lw=lw_std)
        zax.plot (tt_dop, dop["z"], color='r', lw=lw_std)

        gps_x = np.interp (tt_dop, tt_gps, gps["x"])
        gps_y = np.interp (tt_dop, tt_gps, gps["y"])
        gps_z = np.interp (tt_dop, tt_gps, gps["z"])
        delta = np.hypot (gps_x-dop["x"], gps_y-dop["y"], gps_z-dop["z"])
        delta_med = np.median(delta)
        dlt.set_ylim (0, delta_med + args.dd)
        dlt.plot (tt_dop, delta, color='k', lw=lw_std*2)

        # Convert seconds-since-epoch numbers into struct_time objects and then to
        # strings (you can use time.localtime() instead of time.gmtime() to get the
        # time in your local timezone)
        formatter = ticker.FuncFormatter(lambda x, pos: time.strftime('%m-%dT%H:%MZ', time.gmtime(x)))
        time_span = tt_dop[-1] - tt_dop[0]
        if time_span > 86400.0:
            round_interval = 6*3600
        elif time_span > 3600.0:
            round_interval = 300
        else:
            round_interval = 1
        tick_span = round_interval * int ((time_span/5)/round_interval)
        locator = ticker.MultipleLocator(tick_span)
        dlt.xaxis.set_major_formatter(formatter)
        dlt.xaxis.set_major_locator(locator)
        plt.xticks(rotation=30, ha='right')

        xax.set_ylabel ('X [km]')
        yax.set_ylabel ('Y [km]')
        zax.set_ylabel ('Z [km]')
        dlt.set_ylabel ('$\Delta$ [km]')
        dlt.set_xlabel ('UTC time')

        xax.legend(loc='best')

        pp.savefig()

main()
