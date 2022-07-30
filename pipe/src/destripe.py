#! /usr/bin/env python3
import os
import sys
import time
from datetime import datetime
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

def destripe (dst, scd, medval):
    print_message('writing destriping correction to L2 file')
    units = dst['support_data']['fitted_slant_column'].units
    corrected_product = dst.product_type
    try:
        grp = dst['support_data']
        if 'destriping_correction' not in grp.variables:
            dst_des = grp.createVariable('destriping_correction',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_des.long_name = 'destriping correction'.format(corrected_product)
            dst_des.comment = 'xtrack dependent {} slant column destriping correcton'.format(corrected_product)
            dst_des.units = units
        else:
            dst_des = grp['destriping_correction']
            add_history = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')+':destriping correction\n'
            dst.history = '{}{}'.format(dst.history,add_history)
        dst_des[:] = medval
        # leave L2 file's fitted_slant_column unchanged
        scd = scd - medval
    except Exception as e:
        print_message(e)
        print_message('writing destriping correction to L2 file',error=True)

def correct_background (dst, scd, bgrcor):
    print_message('writing background correction to L2 file')
    units = dst['support_data']['fitted_slant_column'].units
    corrected_product = dst.product_type
    try:
        grp = dst['support_data']
        if 'background_correction' not in grp.variables:
            dst_bgr = grp.createVariable('background_correction',np.float32,('xtrack'),fill_value=-1.0e30,zlib=True,complevel=Deflate_Level)
            dst_bgr.long_name = 'background correction'.format(corrected_product)
            dst_bgr.comment = 'xtrack dependent {} slant column background correcton'.format(corrected_product)
            dst_bgr.units = units
        else:
            dst_bgr = grp['background_correction']
            add_history = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')+':background correction\n'
            dst.history = '{}{}'.format(dst.history,add_history)
        dst_bgr[:] = bgrcor
        # leave L2 file's fitted_slant_column unchanged
        scd = scd + bgrcor
    except Exception as e:
        print_message(e)
        print_message('writing background correction to L2 file',error=True)

def apply_corrections (control, corrfile, input_files):

    # Logical to update and write results to L2 files
    yn_L2_des = str_to_bool(control['yn_L2_write_destriping'])
    yn_L2_bgr = str_to_bool(control['yn_L2_write_background'])

    with Dataset(corrfile, 'r') as src:
        medval = src.variables['destriping_correction'][:]
        bgrcor = src.variables['background_correction'][:]

    for fp in input_files:
        try:
            with Dataset(fp,'r+') as dst:
                amf = dst['support_data']['amf'][:]
                scd = dst['support_data']['fitted_slant_column'][:]
                if yn_L2_des:
                    destripe (dst, scd, medval)
                if yn_L2_bgr:
                    correct_background (dst, scd, bgrcor)
                # Save corrected VCDs to L2 file
                dst['product']['vertical_column'][:] = scd/amf
        except Exception as e:
            print_message(e)
            print_message('saving correction(s) to file {}'.format(fp),error=True)

def main():
    pipe_dir = os.getenv ("SDPC_PIPE_DIR")
    if pipe_dir is None:
        default_config = None
    else:
        default_config = os.path.join (pipe_dir, "etc/destripe.yml")

    parser = argparse.ArgumentParser(description='Destripe L2 products')
    parser.add_argument('--corrfile', default=None,
                        help="Correction file path")
    parser.add_argument('--config',  default=default_config,
                        help="Configuration file path")
    parser.add_argument('filenames', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        print_message('cannot find configuration file', error=True)

    if not os.path.isfile(args.corrfile):
        print_message('cannot find correction file: {}'.format(args.corrfile), error=True)

    control = yaml.load(open(args.config),Loader=yaml.BaseLoader)
    print_message('loaded {}'.format(args.config))

    apply_corrections (control, args.corrfile, args.filenames)

if __name__ == "__main__":
    main()
