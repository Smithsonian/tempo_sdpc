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

def plot_spectra (pdf, wavelen, model, spec_scaled, weight, title_string):
    fig = plt.figure(figsize=(9, 6.5))
    plt.plot (wavelen, np.ma.masked_where(spec_scaled <= 0, spec_scaled), 'k')
    plt.plot (wavelen, model, 'r', marker='o', fillstyle='none', linestyle=':')
    plt.xlabel(r'$\lambda$ [nm]')
    plt.ylabel(r'Radiance')
    plt.title (title_string)
    pdf.savefig(fig)
    plt.close(fig)

def process_file (filename, mirror_step):
    dset = netCDF4.Dataset(filename, 'r')
    step_grid = dset.variables["mirror_step"][:]
    step_index = step_grid.tolist().index(mirror_step)
    if (step_grid[step_index] != mirror_step):
        raise LookupError ('mirror_step = {} not found'.format(mirror_step))
    wavelen = dset.variables["wavelen"][step_index,:,:]
    model = dset.variables["model"][step_index,:,:]
    spec_scaled = dset.variables["spec_scaled"][step_index,:,:]
    weight = dset.variables["weight"][step_index,:,:]
    num_xtrack = len(dset.dimensions['xtrack'])
    num_wavelen = len(dset.dimensions['wavelen'])
    dset.close()

    filename_pr = os.path.splitext(filename)
    pdf_filename = '{}.pdf'.format(filename_pr[0])
    print(pdf_filename)
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_filename)

    for xtrack in range(0, num_xtrack, 16):
        title_string = "xtrack={} step={}".format(xtrack, mirror_step)
        # fill values have mask=true
        if (sum(weight[xtrack,:].mask) > num_wavelen/2):
            #print('skipping xtrack={}'.format(xtrack))
            continue
        print('plotting xtrack={}'.format(xtrack))
        plot_spectra (pdf, wavelen[xtrack,:], model[xtrack,:], spec_scaled[xtrack,:], weight[xtrack,:], title_string)

    pdf.close()

parser = argparse.ArgumentParser(description='plot wavecal diagnostics')
parser.add_argument('--filename', help="diagnostic file name")
parser.add_argument('--mirror_step', help="mirror step", type=int)
args = parser.parse_args()

process_file (args.filename, args.mirror_step)

