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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rc('text', usetex=True)

def read_image (filename, path, slab):
    nc = netCDF4.Dataset (filename, 'r')
    var = nc[path]
    num_dims = len(var.dimensions)
    if num_dims == 3:
        img = var[slab,:,:]
    elif num_dims == 2:
        img = var[:,:]
    else:
        print ("*** Error: not supported: {} has dimension {}".format(path, num_dims))
        raise
    nc.close()
    return img

def main():
    parser = argparse.ArgumentParser(description='plot image')
    parser.add_argument('--min', help="min pixel value (for color scale)", default=None, type=float)
    parser.add_argument('--max', help="max pixel value (for color scale)", default=None, type=float)
    parser.add_argument('--outfile', help="output filename")
    parser.add_argument('--index', help="image index (in 3D array)", type=int, default=0)
    parser.add_argument('--rot90', help="rotate image for display, ROT90*90 degrees; ROT90>0 is counter-clockwise, ROT90<0 is clockwise", type=int, default=0)
    parser.add_argument('varpath', help="variable path")
    parser.add_argument('filepath', help="file path")
    #parser.add_argument('step', help="mirror_step index")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    img = read_image (args.filepath, args.varpath, args.index)

    if not args.rot90 == 0:
        img = np.rot90(img, k=args.rot90, axes=(0,1))

    plt_filename = args.outfile
    pdf = PdfPages(plt_filename)

    fig, ax = plt.subplots (1, 1)
    # fig.set_size_inches([7,10])
    plt.subplots_adjust (hspace=0, wspace=0)
    fig = plt.figure(1)
    the_fontsize = 15

    if not args.min is None:
        vmin = args.min
    else:
        vmin = np.quantile(img, 0.025)

    if not args.max is None:
        vmax = args.max
    else:
        vmax = np.quantile(img, 0.975)

    ax.tick_params(labelsize=the_fontsize)
    im = ax.imshow(img, vmin=vmin, vmax=vmax, origin='lower')

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar (im, cax=cax)
    cbar.ax.tick_params(labelsize=the_fontsize)
    cbar.ax.yaxis.offsetText.set_fontsize(the_fontsize)  # Seriously? Good lord.

    title_var = "\\verb|%s[%d]|" % (args.varpath, args.index)

    title_file = "\\verb|%s|" % (os.path.basename(args.filepath))

    #fig.suptitle (title_file, fontsize=the_fontsize)
    ax.set_title (title_file + "\n" + title_var, fontsize=the_fontsize)

    plt.tight_layout()

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    pdf.close()

if __name__ == '__main__':
    main()
