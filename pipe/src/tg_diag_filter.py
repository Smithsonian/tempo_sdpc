#! /usr/bin/env python3
# author: JCH & HQW (originally named collect_o2o2_diagnostics.py)
# collecting solshi & radshi for CLDO4
# to call: 
# conda activate earth
# python tg_diag_filter.py l2fnm --diagfile=fnmdiag --outfile=fnmout
#
# updates:
# hqw removes dependence on log file (202506)
#      baseline & NRT has different log output from sdpc
#      NRT solar fit is saved and used later, log file has no solar fit info
#      baseline still does solar fit every granule and contains solar fit info 

# import libraries
import re
import sys
import os
import argparse
import numpy as np
from netCDF4 import Dataset

# constants
Fill_Value = np.float32(-1.e30)

#===============================================
# functions
# 
#-----------------------------------------------
# parse l2 log file to gather solcal information
# no longer used (202506)
def parse_solcal_log (logfile, nx):
    """
    The goal is to parse lines that contain a subtring that looks like this:
    SOLAR FIT          #   6: hw 1/e =  3.055E-01; e_asy =  0.000E+00; k =  3.526E+00; shift =  2.161E-02; squeeze
    Specifically, we want the values of these variables:
        hw, asy, k, shift
    along with the integer index preceding the colon
    """

    values = {"hw":np.full(nx, Fill_Value), \
             "asy":np.full(nx, Fill_Value), \
               "k":np.full(nx, Fill_Value), \
           "shift":np.full(nx, Fill_Value)}

    with open(logfile,"r") as f:
        for line in f:
            if "SOLAR FIT" not in line:
                continue
            # cut out the line section we want, and split at the ':'
            a = line.index('#')
            b = line.index('squeeze')
            halves = line[a:b-1].split(':')
            # parse xtrack [1:nx] from the first piece,
            # subtracting 1 to get a zero-based index
            xtrack = int(halves[0].strip('#')) - 1
            # strip whitespace from the second piece:
            s = re.sub('[\s+]', '', halves[1])
            # remove extraneous characters, leaving a semicolon-delimited
            # string that contains variable=float pairs
            s = s.replace('1/e','')
            s = s.replace('e_','')
            # split this string on semicolons, extract the float values,
            # and store in the corresponding array.
            fields = s.split(';')
            for v in fields:
                tok = v.split('=')
                if len(tok[0]) > 0:
                    values[tok[0]][xtrack] = float(tok[1])
    return values

#--------------------------------------------
# read fit parameter from diag file
def read_named_fit_parameter (ncfile, name):
    with Dataset (ncfile, 'r') as nc:
        names = nc['fit_parameter_names'][:]

    name_index = None
    i = 0
    for n in names:
        if n.startswith(name):
            name_index = i
            break;
        i += 1
    if name_index is None:
        return None

    with Dataset (ncfile, 'r') as nc:
        values = nc['fit_parameter'][:,:,name_index]

    return values

#-------------------------------------------
# read variable from diag file
def read_l2diag_var(ncfile, varnm):
    with Dataset(ncfile, 'r') as nc:
        values = nc[varnm][:]

    return values

#-------------------------------------------
# read variable from l2 file
def read_l2_var (ncfile,groupname,varnm):
    with Dataset(ncfile, 'r') as nc:
        values = nc[groupname][varnm][:]
    return values

#-------------------------------------------
# write variables
# no longer used (202506)
def write_params (outfile, log_params, radshi, radconvfl):
    # par_map gives the correspondence, key:value, between
    #       key = the output file variable name,
    # and value = the param names parsed from the log file:
    par_map = {'solshi':'shift', 'solhw1e':'hw', 'solasy':'asy', 'solk':'k'}

    (nt,nx) = radshi.shape

    with Dataset (outfile, 'w', clobber=True) as dst:
        dx = dst.createDimension ('xtrack', nx)
        dm = dst.createDimension ('mirror_step', nt)

        v_rshi = dst.createVariable ('radshi', np.float32, ('mirror_step','xtrack'), fill_value=Fill_Value)
        v_rshi[:,:] = radshi[:,:]

        for key, value in par_map.items():
            v = dst.createVariable (key, np.float32, ('xtrack'), fill_value=Fill_Value)
            v[:] = log_params[value]

        v_mdqf = dst.createVariable ('mdqfl', np.int16, ('mirror_step','xtrack'), fill_value=-999)
        v_mdqf[:,:] = radconvfl[:,:]

#-------------------------------------------
# write variables
def write_params2 (outfile, solshi, radshi, radconvfl):

    (nt,nx) = radshi.shape

    with Dataset (outfile, 'w', clobber=True) as dst:
        dx = dst.createDimension ('xtrack', nx)
        dm = dst.createDimension ('mirror_step', nt)

        v_rshi = dst.createVariable ('radshi', np.float32, ('mirror_step','xtrack'), fill_value=Fill_Value)
        v_rshi[:,:] = radshi[:,:]

        v = dst.createVariable ('solshi', np.float32, ('xtrack'), fill_value=Fill_Value)
        v[:] = solshi[:]

        v_mdqf = dst.createVariable ('mdqfl', np.int16, ('mirror_step','xtrack'), fill_value=-999)
        v_mdqf[:,:] = radconvfl[:,:]

#===========================================
# main program
def main():
    parser = argparse.ArgumentParser(description='Collect O2O2 diagnostics')
  #  parser.add_argument('--log', default="log_O2O2.txt",
  #                      help="O2O2 log file name")
    parser.add_argument('--outfile', default="diaglog.nc",
                        help="Output netcdf4 file")
    parser.add_argument('--diagfile', default="O2O2_diag.nc",
                        help="Netcdf4 O2O2 diagnostic file")
    parser.add_argument('filename', 
           help = "L2 O2O2 file after fitting before cldo4")

    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()
    
   # log_file = args.log
   # print('logfile:',log_file)
    diag_file = args.diagfile
    print('diagfile:',diag_file)
    outfile = args.outfile
    print('outfile:',outfile)
    l2_file = args.filename
    print('l2file:',l2_file)

    radshi = read_named_fit_parameter (diag_file, 'shi')
    (nt,nx) = radshi.shape
    solshi = read_l2diag_var(diag_file, 'solcal_shift')
    #log_params = parse_solcal_log (log_file, nx)

    # right after fitting, main_data_quality_flag is in product group
    mdqfl = read_l2_var(l2_file, 'product', 'main_data_quality_flag')

    if radshi is not None:
        #write_params (outfile, log_params, radshi, mdqfl)
        write_params2 (outfile, solshi, radshi, mdqfl)

#===========================================
if __name__ == '__main__':
    main()
