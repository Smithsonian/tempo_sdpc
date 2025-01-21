#! /usr/bin/env python3

import re
import sys
import os
import argparse
import numpy as np
from netCDF4 import Dataset

Fill_Value = np.float32(-1.e30)

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
            s = re.sub(r'[\s+]', '', halves[1])
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
        print("*** Error: variable {} is not in {}:fit_parameter_names".format(name, os.path.basename(ncfile)), flush=True)
        return None

    with Dataset (ncfile, 'r') as nc:
        values = nc['fit_parameter'][:,:,name_index]
    return values

def write_params (outfile, params, radshi):
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
            v[:] = params[value]

def main():
    default_outfile="diaglog.nc"
    parser = argparse.ArgumentParser(description='Collect TG code wavecal diagnostics')
    parser.add_argument('--output', default=default_outfile,
                        help="Netcdf4 output file [default={}]".format(default_outfile))
    parser.add_argument('logfile', help="log file name")
    parser.add_argument('diagfile', help="Netcdf4 diagnostic file")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    log_file = args.logfile
    diag_file = args.diagfile
    outfile = args.output

    radshi = read_named_fit_parameter (diag_file, 'shi')
    if radshi is None:
        sys.exit(0)

    (nt,nx) = radshi.shape
    params = parse_solcal_log (log_file, nx)
    write_params (outfile, params, radshi)

if __name__ == '__main__':
    main()
