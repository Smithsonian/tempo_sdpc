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

plt.rc('text', usetex=True)

plt.rc('lines', linewidth=0.5)

class SMC_Series (object):
    def __init__(self, t=None, proc_meas_x=None, proc_meas_y=None):
        self.time = t
        self.proc_meas_x = proc_meas_x
        self.proc_meas_y = proc_meas_y

    def append(self, s):
        if self.time is None:
            self.time = s.time
            self.proc_meas_x = s.proc_meas_x
            self.proc_meas_y = s.proc_meas_y
        else:
            self.time = np.append (self.time, s.time)
            self.proc_meas_x = np.append (self.proc_meas_x, s.proc_meas_x)
            self.proc_meas_y = np.append (self.proc_meas_y, s.proc_meas_y)

def read_smc_series (filename):
    nc = netCDF4.Dataset (filename, 'r')
    t = nc.variables['time'][:]
    proc_meas_x = nc.variables['proc_meas_x'][:]
    proc_meas_y = nc.variables['proc_meas_y'][:]
    s = SMC_Series (t, proc_meas_x, proc_meas_y)
    nc.close()
    return s

def main():
    parser = argparse.ArgumentParser(description='plot SMC')
    parser.add_argument('--tmin', help="start time [sec]", default=None, type=float)
    parser.add_argument('--tmax', help="end time [sec]", default=None, type=float)
    parser.add_argument('--outfile', help="output filename")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    #parser.add_argument('step', help="mirror_step index")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    smc = SMC_Series()
    for f in args.filenames:
        smc.append(read_smc_series (f))

    t0 = smc.time[0]
    t_last = smc.time[-1]

    if args.tmin:
        tmin = args.tmin
    else:
        tmin = t0

    if tmin < t0:
        print ('WARNING: tmin precedes first data point: t0 = {}'.format(t0))

    if args.tmax:
        tmax = args.tmax
    else:
        tmax = t_last

    if tmax > t_last:
        print ('WARNING: tmax lies beyond the last data point: t_last = {}'.format(t_last))

    plt_filename = args.outfile
    pdf = PdfPages(plt_filename)

    fig, axs = plt.subplots (2, 1, sharex='row')
    plt.subplots_adjust (hspace=0, wspace=0)
    (ax1, ax2) = axs
    fig = plt.figure(1)

    t = smc.time-t0

    ax1.set_xlim (tmin-t0, tmax-t0)
    ax1.plot(t, smc.proc_meas_x)
    ax1.set_ylabel ('proc\_meas\_x [$\mu$rad]')
    ax1.set_title ('SMC, t=[{},{}]'.format(tmin, tmax))

    ax2.set_xlim (tmin-t0, tmax-t0)
    ax2.plot(t, smc.proc_meas_y)
    ax2.set_ylabel ('proc\_meas\_y [$\mu$rad]')
    ax2.set_xlabel ('elapsed time [sec]')

    pdf.savefig(fig)
    plt.close(fig)

    pdf.close()

if __name__ == '__main__':
    main()
