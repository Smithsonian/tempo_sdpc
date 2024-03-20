#! /usr/bin/env python3

import sys
import re
import argparse
from netCDF4 import Dataset as NetCDFFile

def replace_bounding_polygon (name, values, text):
    # find the object
    beg = re.search ('OBJECT[ ]*=[ ]*{}'.format(name), text)
    if (beg == None):
        return text
    b = beg.end()
    end = re.search ('[ \n]*END_OBJECT[ ]*=[ ]*{}'.format(name), text)
    e = end.start()
    # update the number of values
    num_val = re.search ('NUM_VAL[ ]*=', text[b:e])
    num_val_end = re.search ('[ \n]*CLASS[ ]*=', text[b:e])
    nvb = b + num_val.end()
    nve = b + num_val_end.start()
    text = text[:nvb] + ' ' + str(len(values)) + text[nve:]
    # update the values
    value = re.search ('VALUE[ ]*=', text[b:e])
    b += value.end()
    #s = " " + " ".join(text[b:e].split())
    new_s = " ({})".format(", ".join(["{}".format(v) for v in values]))
    return text[:b] + new_s + text[e:]

def update_metfile (bounding_polygon, metfile):

    # Parse bounding polygon string
    pt_strings = bounding_polygon[len("POLYGON"):].lstrip("(").rstrip(")").split(',')
    lats = [float(pt.split(' ')[0]) for pt in pt_strings]
    lons = [float(pt.split(' ')[1]) for pt in pt_strings]
    seqno = [i+1 for i in range(len(lons))]

    # Replace bounding polygon fields in .met file
    with open (metfile, 'r') as fp:
        text = fp.read()

    bdry_dict = {'GRINGPOINTLATITUDE':lats, 'GRINGPOINTLONGITUDE':lons, 'GRINGPOINTSEQUENCENO':seqno}
    for name,value in bdry_dict.items():
        text = replace_bounding_polygon (name, value, text)

    with open (metfile, 'w') as fp:
        fp.write(text)

def update_ncfile (geospatial_bounds, ncfile):
    with NetCDFFile (ncfile, "r+") as nc:
        nc.setncattr ("geospatial_bounds", geospatial_bounds)

def main ():
    parser = argparse.ArgumentParser(description='Propagate netcdf file geospatial_bounds attribute to .met or .nc file')
    parser.add_argument('--src', help="Source netcdf file name")
    parser.add_argument('filename', help="Destination file name", nargs='*', default=None)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    source_ncfile = args.src

    # Read bounding polygon string from netcdf file
    with NetCDFFile(source_ncfile, "r") as nc:
        geospatial_bounds = nc.getncattr ("geospatial_bounds")

    for dest in args.filename:
        if dest.endswith (".met"):
            update_metfile (geospatial_bounds, dest)
        elif dest.endswith (".nc"):
            update_ncfile (geospatial_bounds, dest)
        else:
            print ("Skipping unsupported file type: {}".format(dest))
            continue
        print ("Updated: {}".format(dest))

if __name__ == '__main__':
    main()
