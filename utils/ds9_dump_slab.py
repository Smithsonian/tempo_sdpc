#! /usr/bin/env python3

import os, sys
import numpy as np
from netCDF4 import Dataset
import argparse
import subprocess

def read_slab (filepath, varpath, index):
    with Dataset (filepath, "r") as nc:
        if "/" in varpath:
            grp = nc[os.path.dirname(varpath)]
            name = os.path.basename(varpath)
        else:
            grp = nc
            name = varpath.strip('/')
        v = grp[name]
        if v.ndim == 2:
            array = v[:,:]
        elif v.ndim == 3:
            array = v[index,:,:]
        else:
            print ('*** Error: unsupport array shape (ndim={})'.format(v.ndim))
            array = None

    if np.ma.isMaskedArray (array):
        array = np.ma.getdata(array)

    return array

def main():
    parser = argparse.ArgumentParser(description='Dump array slab to be read by ds9')
    parser.add_argument('--index', default=0,
                        help="Slab index")
    parser.add_argument('--var', default=None,
                        help="NetCDF4 variable path")
    parser.add_argument('--outfile', metavar='PATH', default="arr.dat",
                        help="Output file path")
    parser.add_argument('--display', action='store_true',
                        help="Display image in ds9")
    parser.add_argument('filename', help="NetCDF4 filename")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    infile = args.filename
    outfile = args.outfile
    index = args.index

    ary = read_slab (infile, args.var, args.index)
    if ary is None:
        sys.exit(1)

    ary.astype(ary.dtype, order='C', casting='no').tofile(outfile)

    """
    BITPIX    Numpy Data Type
    8         numpy.uint8 (note it is UNsigned integer)
    16        numpy.int16
    32        numpy.int32
    64        numpy.int64
    -32       numpy.float32
    -64       numpy.float64
    """
    bitpix_table = {'uint8':8, 'int16':16, 'int32':32, 'int64':64, 'float32':-32, 'float64':-64}
    bitpix = bitpix_table[ary.dtype.name]

    shape = ary.shape
    ydim = shape[0]
    xdim = shape[1]

    argv = ['ds9', '-array',
            '{outfile}[xdim={xdim},ydim={ydim},bitpix={bitpix},arch=little,skip=0]'.format(**locals()),
            '-zoom', 'to', 'fit']
    print(' '.join (argv), end='\n')

    if args.display:
        p = subprocess.run (argv)

main()

