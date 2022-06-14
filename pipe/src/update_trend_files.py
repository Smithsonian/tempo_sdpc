#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import numpy as np
import argparse
import pathlib
from netCDF4 import Dataset
from shutil import copyfile

import sqlite3
import signal
from threading import Event

Subgroup_Record_Var = 'mean_time'
Xtrack_Sample_Offset = 25
Xtrack_Sample_Interval = 500
Deflate_Level = 1

Prefix = "trend:"

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stderr, **kwargs)

def logprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stdout, **kwargs)

def trend_file_for_product (infile):
    basename = os.path.basename(infile)
    name_tok = basename.split('_')
    level = name_tok[2]
    nc_prod = Dataset (infile, 'r')
    product_type = nc_prod.product_type
    tok = nc_prod.time_coverage_start.split('-')
    year = tok[0]
    month = tok[1]
    trend_file_basename = 'TEMPO_%s_%s_%s%s.nc' % (product_type, level, year, month)
    nc_prod.close()
    trend_dir = os.path.join ("trend", level, product_type, year)
    return trend_file_basename, product_type, trend_dir

def ensure_dest_dimension_exists (src_dim, dst):
    if src_dim.name in dst.dimensions:
        return
    if src_dim.isunlimited():
        dim_size = None
    else:
        dim_size = len(src_dim)
    dst.createDimension (src_dim.name, src_dim.datatype, dim_size)

def ensure_dest_variable_exists (name, src, dst):
    if name in dst.variables:
        return
    var = src[name]
    for dim_name in var.dimensions:
        ensure_dest_dimension_exists (src.dimensions[dim_name], dst)
    dst.createVariable (name, var.datatype, var.dimensions, zlib=True, complevel=Deflate_Level)
    # copy variable attributes to set any fill value before assigning data
    dst[name].setncatts(var.__dict__)

def append_group_vars_to_existing_file (src, dst, dst_time_var, n):
    for name, var in src.variables.items():
        if dst_time_var.name not in var.dimensions:
            continue
        ensure_dest_variable_exists (name, src, dst)
        num_dims = len(var.shape)
        if num_dims == 1:
            dst[name][n:] = src[name][:]
        elif num_dims == 2:
            dst[name][n:,:] = src[name][:,:]
        elif num_dims == 3:
            dst[name][n:,:,:] = src[name][:,:,:]
        else:
            eprint('*** Error: not implemented: append_file_vars with dimension {}'.format(num_dims))

def walktree(top):
    yield top.groups.values()
    for value in top.groups.values():
        yield from walktree(value)

def ensure_group_exists (name, dst):
    if name not in dst.groups:
        dst_grp = dst.createGroup (name)
    else:
        dst_grp = dst.groups[name]
    return dst_grp

def append_vars_to_existing_file (trend_file, target_file):
    with Dataset(trend_file, 'r') as src, Dataset(target_file, 'a') as dst:
        time_var = dst['time']
        n = time_var.shape[0]
        append_group_vars_to_existing_file (src, dst, time_var, n)
        for children in walktree(src):
            for grp in children:
                dst_grp = ensure_group_exists (grp.name, dst)
                # If a subgroup has a separate variable use that
                # instead of the global 'time'
                if Subgroup_Record_Var in dst_grp.variables:
                    time_var = dst_grp[Subgroup_Record_Var]
                    n = time_var.shape[0]
                append_group_vars_to_existing_file (grp, dst_grp, time_var, n)
    return 0

def append_file_vars (trend_file, target_file):
    if not os.path.isfile(target_file):
        try:
            copyfile (trend_file, target_file)
        except IOError:
            eprint("*** Error copying {} to {}".format (trend_file, target_file))
            return -1
        return 0
    else:
        return append_vars_to_existing_file (trend_file, target_file)

def copy_selected_wavecal_params (src, dst, dst_time_var_name):
    num_xtrack = len(src.dimensions['xtrack'])
    xtrack = [k for k in range(Xtrack_Sample_Offset,num_xtrack,Xtrack_Sample_Interval)]
    num_select = len(xtrack)
    src_var = src.variables['wavecal_params']
    # prep to define the wavecal_params variable
    dst.createDimension ('xtrack', num_select)
    dst.createDimension ('wavecal_par', len(src.dimensions['wavecal_par']))
    dst.createVariable ('xtrack', "i4", ('xtrack',))
    dst['xtrack'][:] = xtrack[:]
    dst.createVariable ('wavecal_params', src_var.datatype, (dst_time_var_name, 'xtrack', 'wavecal_par',), zlib=True, complevel=Deflate_Level)
    dst['wavecal_params'].setncatts(src_var.__dict__)
    selected = src_var[:,:,:]
    dst['wavecal_params'][:,:,:] = selected[:,xtrack,:]

def irr_copy_vars (irr_file, in_trend_file):
    # IRR wavecal params are computed only for the mean spectrum,
    # so we'll need to create a subgroup time variable to hold the mean
    band_names = {'band_540_740_nm', 'band_290_490_nm'}
    with Dataset(irr_file, 'r') as irr, Dataset(in_trend_file, 'r+') as nc:
        src_time_var = irr['time']
        dst_time_var_name = Subgroup_Record_Var
        for band in band_names:
            dst_grp = ensure_group_exists (band, nc)
            dst_grp.createDimension (dst_time_var_name, None)
            dst_grp.createVariable (dst_time_var_name, src_time_var.datatype, (dst_time_var_name,), zlib=True, complevel=Deflate_Level)
            dst_grp[dst_time_var_name].setncatts(src_time_var.__dict__)
            dst_grp[dst_time_var_name][0] = src_time_var[0]
            copy_selected_wavecal_params (irr.groups[band], dst_grp, dst_time_var_name)

def rad_copy_inr_mirror_xy (rad_file, in_trend_file):
    inr_file = rad_file.replace ('TEMPO_RAD_L1', 'TEMPO_INR_L1', 1)
    if not os.path.isfile(inr_file):
        return
    with Dataset(inr_file, 'r') as inr:
        x_var = inr['/SMA/EW']
        y_var = inr['/SMA/NS']
        mirror_x = x_var[:]
        mirror_y = y_var[:]
        datatype = x_var.datatype
    with Dataset(in_trend_file, 'r+') as nc:
        x_var = nc.createVariable ('mirror_x', datatype, ('time',), zlib=True, complevel=Deflate_Level)
        y_var = nc.createVariable ('mirror_y', datatype, ('time',), zlib=True, complevel=Deflate_Level)
        x_var[:] = mirror_x[:]
        y_var[:] = mirror_y[:]

def rad_copy_vars (rad_file, in_trend_file):
    rad_copy_inr_mirror_xy (rad_file, in_trend_file)
    band_names = {'band_540_740_nm', 'band_290_490_nm'}
    with Dataset(rad_file, 'r') as rad, Dataset(in_trend_file, 'r+') as nc:
        for band in band_names:
            dst_grp = ensure_group_exists (band, nc)
            copy_selected_wavecal_params (rad.groups[band], dst_grp, 'time')

def make_samedir_backup (path):
    try:
        copyfile (path, path+'.ORIG')
        return 0
    except IOError:
        eprint ("*** Error making file backup: {}".format(path))
        return -1

def drk_append_vars (drk_file, target_trend_file):
    in_trend_file = os.path.join (os.path.dirname(drk_file), 'trend_params.nc')
    return append_file_vars (in_trend_file, target_trend_file)

def irr_append_vars (irr_file, target_trend_file):
    in_trend_file = os.path.join (os.path.dirname(irr_file), 'trend_params.nc')
    make_samedir_backup (in_trend_file)
    irr_copy_vars (irr_file, in_trend_file)
    return append_file_vars (in_trend_file, target_trend_file)

def rad_append_vars (rad_file, target_trend_file):
    in_trend_file = os.path.join (os.path.dirname(rad_file), 'trend_params.nc')
    make_samedir_backup (in_trend_file)
    rad_copy_vars (rad_file, in_trend_file)
    return append_file_vars (in_trend_file, target_trend_file)

Method_Dict = {
'DRK':drk_append_vars,
'IRR':irr_append_vars,
'RAD':rad_append_vars
}

def update_trend_file (name, *args, **kwargs):
    if name in Method_Dict:
        return Method_Dict[name](*args, **kwargs)
    else:
        eprint('*** Error: unsupported product type: {}'.format(name))
        return -1

def collect_trend_params (path, args_dir, arch_dir):
    trend_file_basename, product_type, trend_dir = trend_file_for_product (path)
    if args_dir == None:
        trend_dir = os.path.join (arch_dir, trend_dir)
    else:
        trend_dir = args_dir
    pathlib.Path(trend_dir).mkdir(parents=True, exist_ok=True)
    trend_file = os.path.join (trend_dir, trend_file_basename)
    status = update_trend_file (product_type, path, trend_file)
    if status != 0:
        eprint ('*** Error processing file: {}'.format(path))
    else:
        logprint ('{} update from: {}'.format(product_type, os.path.basename(path)), flush=True)
    return status

def process_file_list (path_list, args_dir, arch_dir):
    processed_files = []
    if len(path_list) == 0:
        return processed_files
    for path in path_list:
        status = collect_trend_params (path, args_dir, arch_dir)
        if status == 0:
            processed_files = processed_files + [path]
    return processed_files

def table_files_matching_trend_status (cur, table_name, trend_status):
    # verify table existence before attempting a query
    cur.execute ("select count(name) from sqlite_master where type='table' and name='{}'".format(table_name))
    result = cur.fetchone()
    if result[0] == 0:
        return []
    cur.execute ("select path from {} where trend_status == {} order by path".format(table_name, trend_status))
    nc_paths = [item for t in cur.fetchall() for item in t]
    return nc_paths

def connect_database (db_file_path):
    conn = sqlite3.connect (db_file_path)
    conn.execute ("pragma foreign_keys=on")
    return conn

def unprocessed_files (db_file_path, table_names):
    paths = []
    with connect_database (db_file_path) as conn:
        for tbl in table_names:
            paths = paths + table_files_matching_trend_status (conn.cursor(), tbl, 0)
    return paths

def table_name_for_file (filename):
    tok = filename.split('_')
    return '{}_{}'.format(tok[1], tok[2])

def mark_as_processed (db_file_path, path_list):
    if len(path_list) == 0:
        return
    with connect_database (db_file_path) as conn:
        for path in path_list:
            basename = os.path.basename (path)
            table_name = table_name_for_file (basename)
            cur = conn.cursor()
            cur.execute ("update {} set trend_status=1 where filename=\"{}\"".format(table_name, basename))

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGTERM, self.handler)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

def run_as_service (args_dir, arch_dir):
    db_file_path = os.getenv ("SDPC_ARCHIVE_DBFILE")
    if db_file_path == None:
        eprint ('*** Error: SDPC_ARCHIVE_DBFILE is not set')
        sys.exit(1)

    table_names = ["DRK_L1", "IRR_L1", "RAD_L1"]

    sig = Signal_Catcher()

    logprint ("Started", flush=True)
    while not sig.caught():
        filenames = unprocessed_files (db_file_path, table_names)
        processed_files = process_file_list (filenames, args_dir, arch_dir)
        mark_as_processed (db_file_path, processed_files)
        sig.wait(60)

    eprint ("Exiting: caught signal = {}".format(sig.signum))
    return 0

def main():
    parser = argparse.ArgumentParser(description='update product trend files')
    parser.add_argument('--dir', help="trend file directory", default=None)
    parser.add_argument('--service', help="run as a service", action='store_true')
    parser.add_argument('filename', help="netCDF4 file name", nargs='?', default=None)
    """
    if len(sys.argv)==1:
         parser.print_usage(sys.stderr)
         sys.exit(0)
    """
    args = parser.parse_args()

    arch_dir = os.getenv ("SDPC_ARCHIVE_DIR")
    if arch_dir == None:
        eprint ('*** Error: SDPC_ARCHIVE_DIR is not set')
        sys.exit(1)

    if args.service:
        status = run_as_service (args.dir, arch_dir)
    else:
        status = collect_trend_params (args.filename, args.dir, arch_dir)

    if (status != 0):
        sys.exit(1)

if __name__ == '__main__':
    main()
