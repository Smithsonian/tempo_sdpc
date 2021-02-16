#! /usr/bin/env python

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
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

def read_refspec_wavecal (filename):
    dset = netCDF4.Dataset (filename, 'r')
    y = dset.variables["Wavelength"][:]
    s_hr = dset.variables["SolarIrradianceHR"][:]
    s_x = dset.variables["SolarIrradiance"][:,:]
    dset.close()
    return dict([('waves', y), ('irr_hr', s_hr), ('irr_x', s_x)]);

def read_refspec_tracegas (filename):
    (y, s) = np.loadtxt (filename, skiprows=17, dtype='double', unpack=True)
    return dict([('waves', y), ('irr', s)])

def cheb_grid (i, p, i0, n):
    a = i0
    b = i0 + n - 1
    x = (2 * i - a - b) / (b - a)
    y = np.polynomial.chebyshev.chebval (x, p)
    return y;

def read_wavelengths (grp, num, xtrack):
    pvar = grp.variables["wavecal_params"]
    p = pvar[0,xtrack,:]
    start_chan = pvar.getncattr ('start_spectral_channel')
    num_chan = pvar.getncattr ('num_spectral_channels')
    return cheb_grid (np.arange(num, dtype='float'), p, start_chan, num_chan)

def read_irr (filename, group, xtrack):
    dset = netCDF4.Dataset (filename, 'r')
    grp = dset.groups[group]
    s = grp.variables["irradiance"][0,xtrack,:]
    y = read_wavelengths (grp, len(s), xtrack)
    dset.close()
    return dict([('waves', y), ('irr', s)]);

def read_wavecal (filename, num, xtrack):
    dset = netCDF4.Dataset (filename, 'r')
    return read_wavelengths (dset, num, xtrack)

def normalize (x):
    return x / np.max(x)

def plt_overlay (ax, y0, irr0, irr_w0, s):
    #ax.plot (y0, irr0)
    ax.set_xlabel (r'$\lambda$ [nm]')
    ax.plot (y0, irr_w0, linewidth=0.5)
    ax.plot (s['waves'], normalize(s['irr']), 'r', linewidth=0.5)

def main():
    file_refspec_tracegas = '/vex/d2/tempo/sdpc/refdata/trace_gas/data_tempo/TEMPO_solarspec_jqsrt2010.dat'
    file_refspec_wavecal_band1 = '/vex/d2/tempo/sdpc/refdata/instrument/wavcal_refspec_TEMPO_ILS.nc'
    file_refspec_wavecal_band2 = '/vex/d2/tempo/sdpc/refdata/instrument/wavcal_refspec_TEMPO_ILS_Band2.nc'

    refspec_tracegas = read_refspec_tracegas (file_refspec_tracegas)
    refspec_wavecal_band1 = read_refspec_wavecal (file_refspec_wavecal_band1)
    refspec_wavecal_band2 = read_refspec_wavecal (file_refspec_wavecal_band2)

    y0 = refspec_tracegas['waves']
    irr0 = normalize(refspec_tracegas['irr'])
    irr_hr0 = np.interp (y0, refspec_wavecal_band1['waves'], normalize(refspec_wavecal_band1['irr_hr']))

    irr_w0_uv = np.interp (y0, refspec_wavecal_band1['waves'], normalize(refspec_wavecal_band1['irr_x'][0,:]))
    #irr_w0_vis = np.interp (y0, refspec_wavecal_band2['waves'], normalize(refspec_wavecal_band2['irr_x'][0,:]))

    band = 'band_290_490_nm'
    file_irr = '/tmp/TEMPO_irr_L1_V01_20130715T052010Z.nc'

    pdf = matplotlib.backends.backend_pdf.PdfPages('wavecal_check.pdf')

    for xtrack in range(0, 2048, 64):
        irr_uv = read_irr (file_irr, band, xtrack)

        #file_irr_wavecal = '/tmp/TEMPO_irr_L1_V01_20130715T052010Z_wavecal_band_290_490_nm.nc'
        #print ('irr wavelength grids from {}'.format(file_irr_wavecal))
        #irr_uv['waves'] = read_wavecal (file_irr_wavecal, len(irr_uv['waves']), xtrack)

        print ('plotting xtrack={}'.format(xtrack))

        fig, ax = plt.subplots (1, 3, sharey=True)

        ax[0].set_xlim (326, 358)
        ax[0].set_title ('H$_2$CO fit window')
        ax[0].set_ylabel (r'$I(\lambda)/I_\mathrm{max}$')
        plt_overlay (ax[0], y0, irr0, irr_w0_uv, irr_uv)
        ax[0].legend(('reference', 'irradiance'),loc='lower left')

        ax[1].set_xlim (388, 402)
        ax[1].set_title (r'Radiance $\lambda$ cal. window')
        plt_overlay (ax[1], y0, irr0, irr_w0_uv, irr_uv)

        ax[2].set_xlim (424, 452)
        ax[2].set_title (r'NO$_2$ fit window')
        plt_overlay (ax[2], y0, irr0, irr_w0_uv, irr_uv)

        #fig.tight_layout()
        fig.suptitle ('xtrack={}'.format(xtrack))
        pdf.savefig(fig)
        plt.close(fig)

    pdf.close()

main()
