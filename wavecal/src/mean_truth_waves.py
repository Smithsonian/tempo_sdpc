#! /usr/bin/env python3

import netCDF4
import numpy as np

def read_wavelength (dset, grpname):
    grp = dset.groups[grpname]
    y = grp.variables["wavelength"][:,:]
    return y

def print_cheb_params (y2d):
    y = np.mean(y2d, axis=0)
    ymin = y[0]
    ymax = y[-1]
    dy = (ymax - ymin)/ (len(y)-1.0)
    print ("min = {}".format(ymin))
    print ("max = {}".format(ymax))
    print (" dy = {}".format(dy))

def main():
    filename = "/home/houck/test_data/level0/2021feb04/truth_SL6/truth_TEMPO_RAD_L0_V01_20130701T170000Z_S005G01.nc"
    band_uv = "band_290_490_nm"
    band_vis = "band_540_740_nm"

    dset = netCDF4.Dataset (filename, 'r')

    print(band_uv)
    uv = read_wavelength (dset, band_uv)
    print_cheb_params (uv)

    print(band_vis)
    vis = read_wavelength (dset, band_vis)
    print_cheb_params (vis)

    dset.close()

main()

