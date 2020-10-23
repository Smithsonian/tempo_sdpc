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

class IRU_Series (object):
    def __init__(self, t=None, gyro_output=None):
        self.time = t
        self.gyro_output = gyro_output

    def append(self, s):
        if self.time is None:
            self.time = s.time
            self.gyro_output = s.gyro_output
        else:
            self.time = np.append (self.time, s.time)
            self.gyro_output = np.append (self.gyro_output, s.gyro_output, axis=0)

def read_iru_series (filename):
    nc = netCDF4.Dataset (filename, 'r')
    t = nc.variables['time'][:]
    gyro_output = nc.variables['gyro_output'][:,:]
    iru = IRU_Series (t, gyro_output/1.e6)
    nc.close()
    return iru

def main():
    parser = argparse.ArgumentParser(description='plot IRU')
    parser.add_argument('--tmin', help="start time [sec]", default=None, type=float)
    parser.add_argument('--tmax', help="end time [sec]", default=None, type=float)
    parser.add_argument('--outfile', help="output filename")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    #parser.add_argument('step', help="mirror_step index")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    iru = IRU_Series()
    for f in args.filenames:
        iru.append(read_iru_series (f))

    t0 = iru.time[0]
    t_last = iru.time[-1]

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

    fig, axs = plt.subplots (4, 1, sharex='row')
    plt.subplots_adjust (hspace=0, wspace=0)
    fig = plt.figure(1)

    axs[0].set_title ('IRU $/10^6$, t=[{},{}]'.format(tmin, tmax))

    t = iru.time - t0

    for i in range(4):
        ax = axs[i]
        ax.set_xlim (tmin-t0, tmax-t0)
        ax.plot(t, iru.gyro_output[:,i])
        ax.set_ylabel ('gyro\_output[{}]'.format(i))

    axs[3].set_xlabel ('elapsed time [sec]')

    pdf.savefig(fig)
    plt.close(fig)

    pdf.close()

if __name__ == '__main__':
    main()
