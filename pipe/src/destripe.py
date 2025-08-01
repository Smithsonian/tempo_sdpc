#! /usr/bin/env python3
import os
import sys
import time
from datetime import datetime, timezone
import argparse
import yaml
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

Deflate_Level=1

def print_message(message, error=False):
    if error:
        print('{}: ERROR {}'.format(time.asctime(),message))
        sys.exit(1)
    else:
        print('{}: {}'.format(time.asctime(),message))

def str_to_bool(s):
    ''' Convert string (s) to boolean
        ARGS:
         s: string
        RETURNS:
         boolean
    '''
    if s in ('True','true','T','t'):
        return True
    else:
        return False

def destripe (dst, scd, stripe_val):
    units = dst['support_data']['fitted_slant_column'].units
    corrected_product = dst.product_type
    try:
        grp = dst['support_data']
        if 'destriping_correction' not in grp.variables:
            dst_des = grp.createVariable('destriping_correction',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_des.long_name = 'destriping correction'.format(corrected_product)
            dst_des.comment = 'across track dependent {} slant column destriping correcton'.format(corrected_product)
            dst_des.units = units
        else:
            dst_des = grp['destriping_correction']
        add_history = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')+':destriping correction\n'
        dst.history = '{}{}'.format(dst.history,add_history)
        dst_des[:] = np.repeat(stripe_val[np.newaxis,:],scd.shape[0],axis=0)
        # Save original fitted slant column density from spectral fit
        if 'fitted_slant_column_uncorrected' not in grp.variables:
            dst_scd_orig = grp.createVariable('fitted_slant_column_uncorrected',np.float32,('mirror_step','xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_scd_orig.long_name = 'fitted slant column before destriping correction'.format(corrected_product)
            dst_scd_orig.comment = 'fitted slant column before destriping correction'.format(corrected_product)
            dst_scd_orig.coordinates = 'time longitude latitude'
            dst_scd_orig.units = units
            dst_scd_orig[:] = scd

        return (scd - stripe_val)

    except Exception as e:
        print_message(e)
        print_message('destriping L2 file',error=True)

def apply_destripe (corrfile, input_files):

    with Dataset(corrfile, 'r') as src:
        stripe_val = src.variables['destriping_correction'][:]

    for fp in input_files:
        try:
            with Dataset(fp,'r+') as dst:
                print_message('destriping L2 file {}'.format(fp))
                grp_prod = dst['product']
                grp_supp = dst['support_data']
                if 'amf_total' in grp_supp.variables:
                    amf = dst['support_data']['amf_total'][:]
                elif 'amf' in grp_supp.variables:
                    amf = dst['support_data']['amf'][:]
                scd = dst['support_data']['fitted_slant_column'][:]
                scd = destripe (dst, scd, stripe_val)
                # Save corrected SCDs and VCDs to L2 file
                dst['support_data']['fitted_slant_column'][:] = scd
                vcd = scd/amf
                if 'vertical_column' in grp_prod.variables:
                    dst['product']['vertical_column'][:] = vcd
                elif 'vertical_column_total' in grp_supp.variables:
                    dst['support_data']['vertical_column_total'][:] = vcd
        except Exception as e:
            print_message(e)
            print_message('destriping L2 file {}'.format(fp),error=True)

def main():
    parser = argparse.ArgumentParser(description='Destripe L2 products')
    parser.add_argument('--corrfile', default=None,
                        help="Correction file path")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if not os.path.isfile(args.corrfile):
        print_message('cannot find correction file: {}'.format(args.corrfile), error=True)

    apply_destripe (args.corrfile, args.filenames)

if __name__ == "__main__":
    main()
