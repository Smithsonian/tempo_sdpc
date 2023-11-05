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
    def __init__(self, time=None, gyro_output=None, gyro_quality_flag=None,
                       bias_time=None, bias=None):
        self.time = time
        self.gyro_output = gyro_output
        self.gyro_quality_flag = gyro_quality_flag
        self.bias_time = bias_time
        self.bias = bias

    def append(self, s):
        if self.time is None:
            self.time = s.time
            self.gyro_output = s.gyro_output
            self.gyro_quality_flag = s.gyro_quality_flag
            self.bias_time = s.bias_time
            self.bias = s.bias
        else:
            self.time = np.append (self.time, s.time)
            self.gyro_output = np.append (self.gyro_output, s.gyro_output, axis=0)
            self.gyro_quality_flag = np.append (self.gyro_quality_flag, s.gyro_quality_flag, axis=0)
            self.bias_time = np.append (self.bias_time, s.bias_time)
            self.bias = np.append (self.bias, s.bias, axis=0)

    def sort(self):
        it = np.argsort (self.time)
        self.time = self.time[it]
        self.gyro_output = self.gyro_output[it,:]
        self.gyro_quality_flag = self.gyro_quality_flag[it,:]
        ibt = np.argsort (self.bias_time)
        self.bias_time = self.bias_time[ibt]
        self.bias = self.bias[ibt,:]

    def filter_gyro_output(self):
        i = np.nonzero (self.gyro_quality_flag != 0)
        self.gyro_output[i] = 0

def read_iru_series (filename):
    nc = netCDF4.Dataset (filename, 'r')
    time = nc.variables['time'][:]
    gyro_output = nc.variables['gyro_output'][:,:]
    gyro_quality_flag = nc.variables['gyro_quality_flag'][:,:]
    bias_time = nc.variables['bias_time'][:]
    bias = nc.variables['bias'][:,:]
    iru = IRU_Series (time, gyro_output, gyro_quality_flag, bias_time, bias)
    nc.close()
    return iru

def autoscale_ylim (ax, x, y):
    xlim = ax.get_xlim()
    i = np.where ((xlim[0] < x) & (x < xlim[1]))[0]
    ymin = y[i].min()
    ymax = y[i].max()
    pad = 0.02 * (ymax - ymin)
    ax.set_ylim (ymin - pad, ymax + pad)

def process_gyro_axis (iru, i):
    time = iru.time
    output = iru.gyro_output[:,i]

    # At power-on, we uploaded these IRU scale factors:
    #IRU_SCALE_FACT SCALE_FACT_0=0.0037876 SCALE_FACT_1=0.0037876 SCALE_FACT_2=0.0037876 SCALE_FACT_3=0.0037876
    scale_factor = np.array ([0.0037876,0.0037876,0.0037876,0.0037876])
    scale = scale_factor[i];

    dt = np.diff (time)
    dg = np.diff (output)

    # gyro_output is the integrated roll rate, so the
    # value wraps periodically.  Filter out the points
    # where this happens by excluding points with a large
    # discontinuous change in gyro_output:
    max_abs_out = np.max(np.abs(output))
    k = np.nonzero (abs(dg) < max_abs_out)
    dt = dt[k]
    dg = dg[k]

    time = time[k]
    rate = scale * dg/dt

    # apply axis calilbration from Maxar:
    #axis_calibration = np.array([0.053523430394493,0.053523430394493,-0.071364573859324,-0.169490862915894])
    #rate += axis_calibration[i]

    # compute integrated rates
    angle = np.cumsum (rate * dt)

    result = {}
    result["axis"] = i
    result["time"] = time
    result["rate"] = rate
    result["angle"] = angle

    return result

def process_gyro_output (iru):
    axes = []
    for i in range(4):
        axes.append(process_gyro_axis (iru, i))
    return axes

def sc_to_gyro_axes_transformation_matrix ():
    m = np.array(
         [[ 0.600164,  0.551662, 0.579200],   # HRG_A
          [ 0.598222, -0.601944, 0.528955],   # HRG_B
          [-0.555617, -0.601944, 0.573545],   # HRG_C
          [-0.553675,  0.551662, 0.623789]]   # HRG_D
         )
    return m

def main():
    parser = argparse.ArgumentParser(description='plot IRU')
    parser.add_argument('--tmin', help="start time [sec]", default=None, type=float)
    parser.add_argument('--tmax', help="end time [sec]", default=None, type=float)
    parser.add_argument('--outfile', help="output filename", default='iru.pdf')
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    #parser.add_argument('step', help="mirror_step index")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    iru = IRU_Series()
    for f in args.filenames:
        if f[0] != '@':
            iru.append(read_iru_series (f))
        else:
            with open(f[1:]) as fn:
                lines = [line.rstrip() for line in fn]
            for fn in lines:
                iru.append(read_iru_series (fn))

    iru.sort()
    iru.filter_gyro_output

    gyro_axes = process_gyro_output (iru)

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
    print('Creating plot file: {}'.format(plt_filename))
    with PdfPages(plt_filename) as pdf:

        fig, axs = plt.subplots (4, 1, sharex='row')
        plt.subplots_adjust (hspace=1.25, wspace=0)

        plt.suptitle ('Gyro rates [$\mu$rad s${}^{-1}$], t0=%f' % (t0))
        mean_rate = np.zeros(4)

        for g in gyro_axes:
            i = g["axis"]
            t = g["time"] - t0
            rate = g["rate"]
            mean_rate[i] = np.mean(rate)
            ax = axs[i]
            ax.set_xlim (tmin-t0, tmax-t0)
            autoscale_ylim (ax, t, rate)
            ax.plot(t, rate)
            ax.set_ylabel ('axis %d' % (i))
            ax.set_title ('mean = %0.4f $\mu$rad s${}^{-1}$' % (mean_rate[i]), loc='right')
            #if i < 3:
                #    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        axs[3].set_xlabel ('elapsed time [sec]')
        pdf.savefig(fig)
        plt.close(fig)

        fig, axs = plt.subplots (4, 1, sharex='row')
        plt.subplots_adjust (hspace=0, wspace=0)

        plt.suptitle ('Integrated rates [rad], t0=%f' % (t0))

        for g in gyro_axes:
            i = g["axis"]
            t = g["time"] - t0
            angle = g["angle"] * 1.e-6
            ax = axs[i]
            ax.set_xlim (tmin-t0, tmax-t0)
            autoscale_ylim (ax, t, angle)
            ax.plot(t, angle)
            ax.set_ylabel ('axis %d' % (i))
            if i < 3:
                ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

        axs[3].set_xlabel ('elapsed time [sec]')
        pdf.savefig(fig)
        plt.close(fig)

    #pdf.close()

    print("mean gyro rates = ", mean_rate, "urad/sec")

    gyro_axes_used_by_spacecraft = [0,1,2]
    sc_to_gyro = sc_to_gyro_axes_transformation_matrix()
    gyro_to_sc = np.linalg.inv (sc_to_gyro [gyro_axes_used_by_spacecraft, :])

    mean_body_axis_rates = np.matmul (gyro_to_sc, mean_rate[gyro_axes_used_by_spacecraft])
    print("mean body axis rates = ", mean_body_axis_rates, "urad/sec")

    mean_bias = np.mean(iru.bias, axis=0)
    print("mean bias = ", mean_bias)

    corr_mean_body_axis_rate = mean_body_axis_rates - mean_bias;
    print ("bias corrected body axis rates = ", corr_mean_body_axis_rate)

    ideal_sidereal_rates = [0.0, -72.9212, 0.0];  # urad/sec
    print ("ideal sidereal rate: ", ideal_sidereal_rates, "urad/sec")

    print("sidereal subtracted (net) body axis rates = ",
          corr_mean_body_axis_rate - ideal_sidereal_rates, "urad/sec");

if __name__ == '__main__':
    main()
