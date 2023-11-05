#! /usr/bin/env python

import sys
import math
import glob
import argparse

import netCDF4
import numpy as np
#from numpy import ma

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def filter_bad(arr):
    bad = ~np.isfinite(arr)
    num_bad = np.count_nonzero (bad)
    if num_bad > 0:
        num_nan = np.count_nonzero(np.isnan(arr))
        num_inf = np.count_nonzero(np.isinf(arr))
        print('array has {} bad values ({} nan, {} inf)'.format(num_bad, num_nan, num_inf))
    arr[bad] = 1.e-10  # Kind of a hack
    return arr

def rescale(arr):
    arr_min = arr.min()
    arr_max = arr.max()
    if arr_max > arr_min:
        arr_rescaled = (arr - arr_min) / (arr_max - arr_min)
    else:
        arr_rescaled = arr
    return arr_rescaled

def get_image (filename):
    dset = netCDF4.Dataset (filename, 'r')

    grp = dset.groups['band_540_740_nm']
    vis3d = grp.variables['radiance'][:,:,:]
    r_img = filter_bad(vis3d[:,:,0:512].sum(axis=2))
    g_img = filter_bad(vis3d[:,:,513:].sum(axis=2))

    grp = dset.groups['band_290_490_nm']
    uv3d = grp.variables['radiance'][:,:,:]
    b_img = filter_bad(uv3d.sum(axis=2))

    dset.close()

    arr = np.stack ([r_img, g_img, b_img], axis=2)
    return arr

def main():
    parser = argparse.ArgumentParser(description='plot radiances')
    parser.add_argument('--outfile', default='rad.png', help="output filename")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    outfile = args.outfile
    file_list = args.filenames

    num_files = len(file_list)

    file_list.reverse()

    img = get_image (file_list[0])
    print(file_list[0])
    for i in range(1,num_files):
        print(file_list[i])
        img_i = get_image (file_list[i])
        #print(img.shape)
        img = np.vstack ((img_i, img))
        #print(img.shape)

    img = np.rot90(img,1,(1,0))
    img = rescale(img)

    fig = plt.figure(dpi=300, figsize=(3,6))
    plt.tight_layout()
    plt.tick_params(labelsize=5)

    im = plt.imshow(np.sqrt(img), interpolation='none')
    #plt.colorbar(im, orientation='horizontal', labelsize=5)
    plt.savefig (outfile, bbox_inches='tight')

if __name__ == '__main__':
    main()
