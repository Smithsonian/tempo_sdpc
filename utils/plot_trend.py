#! /usr/bin/env python3

import os
import socket
import sys
import subprocess
import math
import argparse
import netCDF4
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib.ticker import ScalarFormatter

plt.rc('text', usetex=True)

plt.rc('lines', linewidth=0.5)

def autoscale_ylim (ax, x, y):
    xlim = ax.get_xlim()
    i = np.where ((xlim[0] < x) & (x < xlim[1]))[0]
    yshape = y.shape
    if len(yshape) == 1:
        ymin = y[i].min()
        ymax = y[i].max()
    else:
        ymin = y[i,:].min()
        ymax = y[i,:].max()
    if ymax > ymin:
        pad = 0.05 * (ymax - ymin)
    else:
        pad = 0.05 * abs(ymin + ymax)/2.0
    ax.set_ylim (ymin - pad, ymax + pad)
    y_formatter = ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)

def annotate_lines (ax, labels):
    for line, name in zip(ax.lines, labels):
        y = line.get_ydata()[-1]
        ax.annotate(name, xy=(1,y), xytext=(6,0), color=line.get_color(),
            xycoords = ax.get_yaxis_transform(), textcoords="offset points",
            size=12, va="center")

def product_type (filename):
    with netCDF4.Dataset (filename, 'r') as nc:
        type_str = nc.product_type
    return type_str

def read_xtrack (filename, grp):
    with netCDF4.Dataset (filename, 'r') as nc:
        xtrack = nc["/%s/xtrack"%(grp)][:]
    return xtrack

class Trend_Series (object):
    def __init__(self, t=None, var=None, time_var_name="time", varpath=None):
        self.time = t
        self.var = var
        self.time_var_name = time_var_name
        self.varpath = varpath

    def append(self, s):
        if self.time is None:
            self.time = s.time
            self.var = s.var
        else:
            self.time = np.append (self.time, s.time)
            self.var = np.append (self.var, s.var, axis=0)

    def sort(self):
        i = np.argsort(self.time)
        self.time = self.time[i]
        num_dims = len(self.var.shape)
        if num_dims == 1:
            self.var = self.var[i]
        elif num_dims == 2:
            self.var = self.var[i,:]
        elif num_dims == 3:
            self.var = self.var[i,:,:]

def read_trend_series (filename, varpath):
    nc = netCDF4.Dataset (filename, 'r')

    # the first dimension must specify the time variable
    nc_var = nc[varpath]
    var_dims = nc_var.dimensions
    time_var_name = var_dims[0]

    # the time variable might be in the same group
    # or it might be at the top level
    if "/" in varpath:
        grp_path = os.path.dirname(varpath)
        grp_variables = nc[grp_path].variables.keys()
    else:
        grp_path = "/"
        grp_variables = nc.variables.keys()

    if time_var_name in grp_variables:
        timevarpath = os.path.join (grp_path, time_var_name)
    else:
        timevarpath = "time"
    t = nc[timevarpath][:]

    # the variable to plot may have more than 1 dimension:
    num_dims = len(nc_var.shape)
    if num_dims == 1:
        var = nc_var[:]
    elif num_dims == 2:
        var = nc_var[:,:]
    elif num_dims == 3:
        var = nc_var[:,:,:]
    else:
        print("*** Error: unsupported variable shape: num_dims = {}".format(num_dims))
        sys.exit(1)
    trend = Trend_Series (t, var, time_var_name, varpath)
    nc.close()
    return trend

def read_trend (varpath, filenames):
    trend = Trend_Series(varpath=varpath)
    for f in filenames:
        trend.append(read_trend_series (f, varpath))
    trend.sort()
    return trend

def filter_time_range (tmin, tmax, trend):
    t0 = trend.time[0]
    t_last = trend.time[-1]
    if tmin == None:
        tmin = t0
    else:
        tmin = tmin
    if tmin < t0:
        print ('WARNING: tmin precedes first data point: t0 = {}'.format(t0))
    if tmax == None:
        tmax = t_last
    else:
        tmax = tmax
    if tmax > t_last:
        print ('WARNING: tmax lies beyond the last data point: t_last = {}'.format(t_last))
    return tmin, tmax

def plot_eoffsets (pdf, filenames, in_tmin, in_tmax):
    trend = read_trend ("eoffsets", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, trend)
    t = trend.time - tmin
    fig = plt.figure()
    gs = fig.add_gridspec(1,1, hspace=0)
    ax = fig.add_subplot(gs[0])
    ax.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax, t, trend.var)
    ax.plot(t, trend.var)
    ax.set_ylabel ("eoffset")
    ax.set_xlabel ("time [sec]")
    ax.set_title ("electronic offsets")
    annotate_lines (ax, ['Ao', 'Bo', 'Co', 'Do', 'Ae', 'Be', 'Ce', 'De'])
    pdf.savefig(fig)
    plt.close(fig)

def plot_gain (pdf, filenames, in_tmin, in_tmax):
    trend = read_trend ("gain", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, trend)
    t = trend.time - tmin
    fig = plt.figure()
    gs = fig.add_gridspec(1,1, hspace=0)
    ax = fig.add_subplot (gs[0])
    ax.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax, t, trend.var)
    ax.plot(t, trend.var)
    ax.set_ylabel ("gain [electrons/DN]")
    ax.set_xlabel ("time [sec]")
    ax.set_title ("gain")
    annotate_lines (ax, ['Ao', 'Bo', 'Co', 'Do', 'Ae', 'Be', 'Ce', 'De'])
    pdf.savefig(fig)
    plt.close(fig)

def plot_temps (pdf, filenames, in_tmin, in_tmax):
    fpa_temp = read_trend ("fpa_temp", filenames)
    fpe_temp = read_trend ("fpe_temp", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, fpa_temp)
    t_fpa = fpa_temp.time - tmin
    t_fpe = fpe_temp.time - tmin

    fig = plt.figure()
    gs = fig.add_gridspec(2,1, hspace=0)
    ax1 = fig.add_subplot (gs[0])
    ax2 = fig.add_subplot (gs[1], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax1, t_fpa, fpa_temp.var)
    ax1.tick_params(direction="in", bottom=True, top=True)
    ax1.plot(t_fpa, fpa_temp.var)
    ax1.set_ylabel ("FPA [C]")
    ax1.set_title ("Temperature")

    ax2.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax2, t_fpe, fpe_temp.var)
    ax2.plot(t_fpe, fpe_temp.var)
    ax2.set_ylabel ("FPE [C]")
    ax2.set_xlabel ("time [sec]")
    ax2.tick_params(direction="in", top=True)

    pdf.savefig(fig)
    plt.close(fig)

def plot_solar_angles (pdf, filenames, in_tmin, in_tmax):
    theta = read_trend ("solar_theta", filenames)
    phi = read_trend ("solar_phi", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, theta)
    t = theta.time - tmin

    fig = plt.figure()
    gs = fig.add_gridspec(2,1, hspace=0)
    ax1 = fig.add_subplot (gs[0])
    ax2 = fig.add_subplot (gs[1], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax1, t, theta.var)
    ax1.tick_params(direction="in", bottom=True, top=True)
    ax1.plot(t, theta.var)
    ax1.set_ylabel (r"$\theta$ [deg]")
    ax1.set_title ("solar boresight angles")

    ax2.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax2, t, phi.var)
    ax2.plot(t, phi.var)
    ax2.set_ylabel (r"$\phi$ [deg]")
    ax2.set_xlabel ("time [sec]")
    ax2.tick_params(direction="in", top=True)

    pdf.savefig(fig)
    plt.close(fig)

def plot_sdc (pdf, filenames, in_tmin, in_tmax):
    trend = read_trend ("dc_storage_region", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, trend)
    t = trend.time - tmin
    fig = plt.figure()
    gs = fig.add_gridspec(4,1, hspace=0)
    ax0 = fig.add_subplot(gs[0,0])
    ax0.set_title ("storage region dark current [electrons/s]")
    axe = fig.add_subplot(gs[-1,0])
    axe.set_xlabel ("time [s]")
    quad_labels = ['A', 'B', 'C', 'D']
    for i in range(4):
        ax = fig.add_subplot (gs[i,0])
        ax.set_xlim (0.0, tmax-tmin)
        autoscale_ylim (ax, t, trend.var[:,i])
        #ax.plot(t, trend.var[:,i], linestyle="None", marker=',', alpha=0.2
        ax.hist2d(t, trend.var[:,i], (100, 32), cmap=plt.cm.Blues, linewidth=0, rasterized=True)
        if i < 3:
            plt.setp(ax.get_xticklabels(), visible=False)
        #ax.set_ylabel (quad_labels[i], rotation=0)
        ax.annotate (quad_labels[i], (0.05, 0.8), xycoords = "axes fraction")
    pdf.savefig(fig)
    plt.close(fig)

def plot_dark (pdf, filenames, in_tmin, in_tmax):
    mean = read_trend ("dc_mean", filenames)
    stddev = read_trend ("dc_stddev", filenames)
    tmin, tmax = filter_time_range (in_tmin, in_tmax, mean)
    t = mean.time - tmin
    fig = plt.figure()
    gs = fig.add_gridspec(2,1, hspace=0)
    ax1 = fig.add_subplot (gs[0])
    ax2 = fig.add_subplot (gs[1], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax1.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax1, t, mean.var)
    ax1.tick_params(direction="in", bottom=True, top=True)
    ax1.plot(t, mean.var)
    ax1.set_ylabel ("mean DC [electrons/s]")
    ax1.set_title ("dark current")
    annotate_lines (ax1, ['A', 'B', 'C', 'D'])

    ax2.set_xlim (0.0, tmax-tmin)
    autoscale_ylim (ax2, t, stddev.var)
    ax2.plot(t, stddev.var)
    ax2.set_ylabel ("DC std. dev. [electrons/s]")
    ax2.set_xlabel ("time [s]")
    annotate_lines (ax2, ['A', 'B', 'C', 'D'])
    ax2.tick_params(direction="in", top=True)

    pdf.savefig(fig)
    plt.close(fig)

def plot_pqf_bits (pdf, trend, in_tmin, in_tmax):
    tmin, tmax = filter_time_range (in_tmin, in_tmax, trend)
    t = trend.time - tmin

    shape = trend.var.shape
    num_bits = shape[1]
    bits_used = []
    for i in range(num_bits):
        if any(trend.var[:,i] > 0):
            bits_used.append(i)

    num_bits_used = len(bits_used)
    if num_bits_used == 0:
        return

    bit_labels = ["missing data", "bad pixel", "processing error", "transient pixel", "rts pixel",
                  "saturated", "noise underflow", "dark corr. error", "offset corr. error", "smear corr. error",
                  "straylight corr. error", "nonlinear range error", "hot pixel", "cold pixel", "unused1", "unused2"]

    fig = plt.figure()
    gs = fig.add_gridspec(num_bits_used,1, hspace=0)
    for i in range(num_bits_used):
        ax = fig.add_subplot (gs[i,0])
        bit = bits_used[i]
        ax.set_xlim (0.0, tmax-tmin)
        autoscale_ylim (ax, t, trend.var[:,bit])
        ax.plot(t, trend.var[:,bit])
        ax.annotate (bit_labels[bit], (0.75, 0.75), xycoords="axes fraction")
        if i < num_bits_used-1:
            plt.setp(ax.get_xticklabels(), visible=False)

    ax0 = fig.add_subplot(gs[0,0])
    ax0.set_title (r"\verb|%s|"%(os.path.basename(trend.varpath)))
    axe = fig.add_subplot(gs[-1,0])
    axe.set_xlabel ("time [sec]")
    pdf.savefig(fig)
    plt.close(fig)

def plot_wavecal_band (pdf, params, xtrack, in_tmin, in_tmax):
    tmin, tmax = filter_time_range (in_tmin, in_tmax, params)
    if tmin == tmax:
        return
    t = params.time - tmin
    fig = plt.figure()
    gs = fig.add_gridspec(2,1, hspace=0)
    ax1 = fig.add_subplot (gs[0])
    ax2 = fig.add_subplot (gs[1], sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)

    n = len(xtrack)

    ax1.set_xlim (0.0, tmax-tmin)
    ax1.set_title (r"\verb|%s|"%(params.varpath))
    ax1.set_ylabel ("$c_0$ [nm]")
    ax1.tick_params(direction="in", bottom=True, top=True)
    ax2.set_ylabel ("$c_1$ [nm/pixel]")
    ax2.set_xlabel ("time [s]")
    ax2.tick_params(direction="in", top=True)

    plt.setp(ax1.get_xticklabels(), visible=False)

    autoscale_ylim (ax1, t, params.var[:,:,0])
    for i in range(n):
        ax1.plot(t, params.var[:,i,0])

    autoscale_ylim (ax2, t, params.var[:,:,1])
    for i in range(n):
        ax2.plot(t, params.var[:,i,1])

    """
    xtrack_labels = ["%d"%(xtrack[i]) for i in range(n)]
    annotate_lines (ax1, xtrack_labels)
    annotate_lines (ax2, xtrack_labels)
    """

    pdf.savefig(fig)
    plt.close(fig)

def plot_wavecal (pdf, filenames, in_tmin, in_tmax):
    xtrack = read_xtrack (filenames[0], "band_290_490_nm")
    wave_par_uv = read_trend ("/band_290_490_nm/wavecal_params", filenames)
    plot_wavecal_band (pdf, wave_par_uv, xtrack, in_tmin, in_tmax)
    xtrack = read_xtrack (filenames[0], "band_540_740_nm")
    wave_par_vis = read_trend ("/band_540_740_nm/wavecal_params", filenames)
    plot_wavecal_band (pdf, wave_par_vis, xtrack, in_tmin, in_tmax)

def main():
    parser = argparse.ArgumentParser(description='plot trend parameters')
    parser.add_argument('--tmin', help="start time [sec]", default=None, type=float)
    parser.add_argument('--tmax', help="end time [sec]", default=None, type=float)
    parser.add_argument('--outfile', help="output filename")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    plt_filename = args.outfile

    type_str = product_type (args.filenames[0])

    with PdfPages(plt_filename) as pdf:
        plot_eoffsets (pdf, args.filenames, args.tmin, args.tmax)
        plot_gain (pdf, args.filenames, args.tmin, args.tmax)
        plot_temps (pdf, args.filenames, args.tmin, args.tmax)
        plot_sdc (pdf, args.filenames, args.tmin, args.tmax)

        if type_str == "TREND_IRR":
            plot_solar_angles (pdf, args.filenames, args.tmin, args.tmax)

        if type_str == "TREND_DRK":
            plot_dark (pdf, args.filenames, args.tmin, args.tmax)
            pqf_bits = read_trend ("pqf_bits", args.filenames)
            plot_pqf_bits (pdf, pqf_bits, args.tmin, args.tmax)

        if type_str == "TREND_IRR" or type_str == "TREND_RAD":
            pqf_bits_uv = read_trend ("pqf_bits_uv", args.filenames)
            plot_pqf_bits (pdf, pqf_bits_uv, args.tmin, args.tmax)
            pqf_bits_vis = read_trend ("pqf_bits_vis", args.filenames)
            plot_pqf_bits (pdf, pqf_bits_vis, args.tmin, args.tmax)
            plot_wavecal (pdf, args.filenames, args.tmin, args.tmax)

if __name__ == '__main__':
    main()
