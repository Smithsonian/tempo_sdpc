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

class HK_Series (object):
    def __init__(self, t=None, var=None):
        self.time = t
        self.var = var

    def append(self, s):
        if self.time is None:
            self.time = s.time
            self.var = s.var
        else:
            self.time = np.append (self.time, s.time)
            self.var = np.append (self.var, s.var)

def read_hk_series (filename, varpath, timevar):
    nc = netCDF4.Dataset (filename, 'r')
    if timevar == None:
        timevarpath = "{}/time".format(os.path.dirname(varpath))
    else:
        timevarpath = timevar
    t = nc[timevarpath][:]
    var = nc[varpath][:]
    hk = HK_Series (t, var)
    nc.close()
    return hk

def autoscale_ylim (ax, x, y):
    xlim = ax.get_xlim()
    i = np.where ((xlim[0] < x) & (x < xlim[1]))[0]
    ymin = y[i].min()
    ymax = y[i].max()
    pad = 0.02 * (ymax - ymin)
    ax.set_ylim (ymin - pad, ymax + pad)
    y_formatter = ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)

def main():
    parser = argparse.ArgumentParser(description='plot HK')
    parser.add_argument('--tmin', help="start time [sec]", default=None, type=float)
    parser.add_argument('--tmax', help="end time [sec]", default=None, type=float)
    parser.add_argument('--outfile', help="output filename")
    parser.add_argument('--timevar', help="time variable path", default=None)
    parser.add_argument('varpath', help="variable path")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    #parser.add_argument('step', help="mirror_step index")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    hk = HK_Series()
    for f in args.filenames:
        hk.append(read_hk_series (f, args.varpath, args.timevar))

    t0 = hk.time[0]
    t_last = hk.time[-1]

    if args.tmin == None:
        tmin = t0
    else:
        tmin = args.tmin

    if tmin < t0:
        print ('WARNING: tmin precedes first data point: t0 = {}'.format(t0))

    if args.tmax == None:
        tmax = t_last
    else:
        tmax = args.tmax

    if tmax > t_last:
        print ('WARNING: tmax lies beyond the last data point: t_last = {}'.format(t_last))

    plt_filename = args.outfile
    pdf = PdfPages(plt_filename)

    fig, axs = plt.subplots (1, 1)
    plt.subplots_adjust (hspace=0, wspace=0)
    fig = plt.figure(1)

    ax = axs
    
    varpath_label = "\\verb|{}|".format(args.varpath)
    fig.suptitle (varpath_label)
    ax.set_title ('t=[{},{}]'.format(tmin, tmax))

    t = hk.time - t0

    ax.set_xlim (tmin-t0, tmax-t0)
    autoscale_ylim (ax, t, hk.var[:])
    ax.plot(t, hk.var[:])
    ax.set_ylabel (varpath_label)
    ax.set_xlabel ('elapsed time [sec]')

    pdf.savefig(fig)
    plt.close(fig)

    pdf.close()

if __name__ == '__main__':
    main()
