#! /usr/bin/env python3

import os
import sys
import re
import datetime
import argparse
from string import Template

import yaml
from netCDF4 import Dataset as NetCDFFile

Sdpc_Root_Dir = os.getenv ('SDPC_ROOT')
if Sdpc_Root_Dir:
    Sdpc_Metadata_Dir = os.path.join (Sdpc_Root_Dir, 'etc/metadata')

def read_keyword_file (filename):
    if not os.path.isfile(filename):
        return {}
    with open(filename, 'r') as stream:
        try:
            meta = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return meta

def write_netcdf_keys (meta, filename, group=None):
    with NetCDFFile(filename, "r+") as nc:
        #nc.set_ncstring_attrs (True)
        if (group == None):
            grp = nc
        elif group in nc.groups:
            grp = nc.groups[group]
        else:
            grp = nc.createGroup(group)
        grp.setncatts (meta)

def write_netcdf_coremetadata (filename, str):
    with NetCDFFile(filename, "r+") as nc:
        nc.setncattr_string ("coremetadata", str)

def read_odl (granule_met_file):
    with open(granule_met_file, 'r') as stream:
        odl = stream.read()
    return odl

def make_template_string (meta):
    # FIXME: this formats everything as a quoted string, but when we
    # receive a proper ODL template, this won't be needed at all.
    fmt = '''
       OBJECT     = %s
          NUM     = 1
        VALUE     = "$%s"
       END OBJECT = %s
    '''
    template = ''
    for k in meta.keys():
        delete_tbl = k.maketrans({'_':None})
        name = k.upper().translate(delete_tbl)
        template += fmt % (name, k, name)
    return template

def insert_odl (file_odl_str, meta):
    odl_tmpl = make_template_string (meta)
    odl_str = Template(odl_tmpl).substitute(meta)
    p = re.compile('(END_GROUP\s*=\s*ADDITIONALATTRIBUTES)')
    split_fstr = p.split (file_odl_str)
    return ''.join([split_fstr[0], odl_str, split_fstr[1], split_fstr[2]])

def write_ascii (filename, str):
    with open(filename,'w') as stream:
        stream.write(str)

def get_day_of_year (nc):
    tstart_str = nc.getncattr('time_coverage_start')
    ymd_str = tstart_str.split('T')[0]
    dt = datetime.datetime.strptime (ymd_str, '%Y-%m-%d')
    day_of_year = dt.timetuple().tm_yday
    return day_of_year

def get_geospatial_bounds (nc):
    if not 'metadata' in nc.groups:
        return {}
    grp = nc.groups['metadata']
    attr = grp.__dict__
    if not 'polygon_longitudes' in attr:
        return {}
    poly_lon = grp.getncattr('polygon_longitudes')
    poly_lat = grp.getncattr('polygon_latitudes')
    bounds_dict = {
    'longitude_min':poly_lon.min(),
    'longitude_max':poly_lon.max(),
    'latitude_min':poly_lat.min(),
    'latitude_max':poly_lat.max()
    }
    return bounds_dict

def collect_keyword_groups (mission_meta, product_meta):
    # Keywords may be labeled as either "global" or "metadata".
    # When no such label is present, the default is "metadata".

    meta_global = {}
    meta_group = {}

    for d in [mission_meta, product_meta]:
        use_default = True
        if 'global' in d:
            if d['global']:
                meta_global.update(d['global'])
                use_default = False
        if 'metadata' in d:
            if d['metadata']:
                meta_group.update(d['metadata'])
                use_default = False
        if use_default and d:
            meta_group.update(d)

    return meta_global, meta_group

def metadata_file_path (basename):
    if not Sdpc_Metadata_Dir:
        return basename
    return os.path.join (Sdpc_Metadata_Dir, basename)

def process_odl (ncfile, meta):
    met_file = ncfile + ".met"
    if not os.path.isfile(met_file):
        return
    odl = read_odl (met_file)
    if False:
        odl = insert_odl (odl, meta)
        write_ascii (met_file + '.new', odl)
    write_netcdf_coremetadata (ncfile, odl)

def main():
    parser = argparse.ArgumentParser(description='insert fixed metadata keywords')
    parser.add_argument('ncfile', help="netCDF data file name")
    args = parser.parse_args()

    ncfile = args.ncfile

    ncfile_prefix = os.path.basename(ncfile.split('_V')[0])

    # Find metadata keyword files for this product
    product_yaml = metadata_file_path (ncfile_prefix + '.yaml')
    mission_yaml = metadata_file_path ("TEMPO.yaml")

    # Read metadata keyword files:
    mission_meta = read_keyword_file (mission_yaml)
    product_meta = read_keyword_file (product_yaml)
    meta_global, meta_group = collect_keyword_groups (mission_meta, product_meta)

    # Derive metadata keywords from file content:
    with NetCDFFile(ncfile, "r") as nc:
        meta_geo = {'day_of_year': get_day_of_year(nc)}
        bounds = get_geospatial_bounds (nc)
        if bounds:
            meta_geo.update()
    meta_global.update(meta_geo)

    # Write metadata to netcdf4 file
    write_netcdf_keys (meta_global, ncfile)
    write_netcdf_keys (meta_group, ncfile, group='metadata')

    # Merge keyword dictionaries:
    meta = {}
    if meta_global:
        meta.update(meta_global)
    if meta_group:
        meta.update(meta_group)

    process_odl (ncfile, meta)

if __name__ == '__main__':
    main()
