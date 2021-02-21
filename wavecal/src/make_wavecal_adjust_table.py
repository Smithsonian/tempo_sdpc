#! /usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import os

"""
When a linear wavelength grid is defined by Chebyshev series
coefficients, the coefficients c[0] and c[1] correspond to:

  c[0] = the midpoint wavelength
       = \lambda( (k0+k1)/2 )
  c[1] = the half-width of the calibration band
       = 0.5 * (\lambda(k1) - \lambda(k0))

where the calibration band is:
   \lambda(k0) <= \lambda <= \lambda(k1)
and k0, k1 are pixel indices.

"""

class Cheb_Params_Type (object):
    def __init__ (self, num_params, a, b, ia, ib, units=None):
        mid = 0.5 * (a + b)
        half_width = 0.5 * (b - a)
        params = np.zeros(num_params)
        params[0] = mid
        params[1] = half_width
        self.params = params
        self.ia = ia
        self.ib = ib
        self.units = units
        self.num_params = num_params

class Entry_Type (object):
    def __init__ (self, name, c_uv, c_vis):
        self.name = name
        self.coeff = {'band_290_490_nm':c_uv, 'band_540_740_nm':c_vis}

def create_band_group (rootgrp, group_name, x, var_list):
    max_num_params = 5
    grp = rootgrp.createGroup (group_name)
    for v in var_list:
        dimname = 'max_num_' + v.name
        vc = v.coeff[group_name]
        dim = grp.createDimension (dimname, max_num_params)
        var = grp.createVariable (v.name, "f4", (x.name, dimname))
        if (vc.units != None):
            var.setncattr ('comment', 'wavelength[i], pix_min <= i <= pix_max')
            var.setncattr ('units', vc.units)
        var.setncattr ('num_params', vc.num_params)
        var.setncattr ('pix_min', vc.ia)
        var.setncattr ('pix_max', vc.ib)
        var_value = np.zeros([x.size, max_num_params])
        var_value[:,0] = vc.params[0]
        var_value[:,1] = vc.params[1]
        var[:,:] = var_value[:,:]

def create_file (filename, var_list):
    ds = Dataset(filename, "w", format="NETCDF4")
    x_name = 'xtrack'
    num_x = 2048
    x = ds.createDimension (x_name, num_x)
    create_band_group (ds, 'band_290_490_nm', x, var_list)
    create_band_group (ds, 'band_540_740_nm', x, var_list)
    ds.close()

class Grid_Type (object):
    def __init__ (self, y0, delta, num):
        self.y = y0 + delta * np.arange(num)
        self.y0 = y0
        self.delta = delta
        self.num = num

def main():

    # UV, Vis wavelength grids:
    num_waves = 1028  # 1026
    delta_lam = 0.198 # 0.2
    uv0 = 291.8       # 290.0
    vis0 = 537.2      # 540.0

    # number of Chebyshev series coefficients
    num_cheb_terms = 2

    # UV narrow band
    iuv0 = 513   # 508  # 510
    nuv0 = 61

    # Vis narrow band
    ivis0 = 242  # 259 # 245
    nvis0 = 61   # 31

    uv = Grid_Type (uv0, delta_lam, num_waves)
    vis = Grid_Type (vis0, delta_lam, num_waves)

    full = {
    'uv':Cheb_Params_Type(num_cheb_terms, uv.y[0], uv.y[-1], 0, num_waves-1, units='nm'),
    'vis':Cheb_Params_Type(num_cheb_terms, vis.y[0], vis.y[-1], 0, num_waves-1, units='nm')
    }

    # Note that the adjustment Chebyshev series _must_ be defined on the same pixels
    # as the full-band Chebyshev series (must have the same x coordinate).
    adjust = {
    'uv':Cheb_Params_Type(num_cheb_terms, 0.0, 0.0, 0, num_waves-1),
    'vis':Cheb_Params_Type(num_cheb_terms, 0.0, 0.0, 0, num_waves-1)
    }

    # For now, shift will only adjust the constant term
    adjust['uv'].params[[0,1]] = [1.0,0.0];
    adjust['vis'].params[[0,1]] = [1.0,0.0];

    iuv1 = iuv0 + nuv0 - 1
    ivis1 = ivis0 + nvis0 - 1

    narrow = {
    'uv':Cheb_Params_Type(num_cheb_terms, uv.y[iuv0], uv.y[iuv1], iuv0, iuv1, units='nm'),
    'vis':Cheb_Params_Type(num_cheb_terms, vis.y[ivis0], vis.y[ivis1], ivis0, ivis1, units='nm')
    }

    var_list = [];
    var_list.append (Entry_Type ('full_band', full['uv'], full['vis']))
    var_list.append (Entry_Type ('full_band_adjust', adjust['uv'], adjust['vis']))
    var_list.append (Entry_Type ('narrow_band', narrow['uv'], narrow['vis']))

    #filename = '/tmp/wavecal_adjust_20190124.nc'
    #filename = '/tmp/wavecal_adjust_20190727.nc'
    filename = '/tmp/wavecal_adjust_20210215.nc'

    create_file (filename, var_list)

main()
