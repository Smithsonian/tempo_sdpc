#! /usr/bin/env python3

import os
import sys
import re
import argparse
import json

def format_CollectionReference (shortname, version):
    s = '"CollectionReference": {"ShortName": %s, "Version": "V%02d"}' % (shortname, int(version))
    return s

def format_DataGranule (granuleid, productiondatetime):
    s = '"DataGranule": {"Identifiers": [{"Identifier": %s,' % (granuleid) \
      + ' "IdentifierType": "ProducerGranuleId"}, {"Identifier": "RFC1321 MD5 = not yet calculated", "IdentifierType": "LocalVersionId"}],' \
      + '"ProductionDateTime": %s}' % (productiondatetime)
    return s

def format_GranuleUR (granuleid):
    return '"GranuleUR": {}'.format(granuleid)

def format_InputGranules (inputpointer):
    return '"InputGranules": [ %s ]' % (re.sub ('[()]', '', inputpointer))

def format_PGEVersionClass (pgeversion):
    return '"PGEVersionClass": {"PGEVersion": %s}' % (pgeversion)

def format_TemporalExtent (beg_date, beg_time, end_date, end_time):
    beg_date = beg_date.strip('"')
    beg_time = beg_time.strip('"')
    end_date = end_date.strip('"')
    end_time = end_time.strip('"')
    s = '"TemporalExtent": { "RangeDateTime": {' \
      + '"BeginningDateTime": "{}T{}",'.format(beg_date, beg_time) \
      + '"EndingDateTime": "{}T{}"'.format(end_date, end_time) \
      + '}}'
    return s

def format_SpatialExtent (lat_str, lon_str):
    lats = re.sub('[()]', '', lat_str).split(',')
    lons = re.sub('[()]', '', lon_str).split(',')

    pts = []
    for lat,lon in zip(lats,lons):
        pts.append ('{"Latitude":%s, "Longitude": %s}' % (lat,lon))

    s = '"SpatialExtent": {"HorizontalSpatialDomain": {"Geometry": {"GPolygons": [{"Boundary": {"Points": [' \
        + ",".join (pts) \
        + ']}}]}}}'

    return s

MetadataSpecification = '"MetadataSpecification": {"Name": "UMM-G", "URL": "https://cdn.earthdata.nasa.gov/umm/granule/v1.6.5", "Version": "1.6.5"}'

def convert_odl_file (metfile, jsonfile=None):

    # Default output filename has .met extension replaced with .cmr.json
    if jsonfile is None:
        jsonfile = os.path.splitext (metfile)[0] + '.cmr.json'

    # Read the input ODL file
    with open (metfile, "r") as fp:
        odl = fp.read()

    # Parse the ODL objects to extract metadata tokens
    object_regex = "OBJECT\s+=\s+(?P<object_name>\w+)\s+NUM_VAL\s+=\s+(?P<num_val>\d+)(?:|\s+CLASS\s+=\s+..*)\s+VALUE\s+=\s+(?P<value>..*)\s+END_OBJECT"
    objc = re.compile (object_regex)
    obj_iter = objc.finditer (odl)

    tokens = {}
    for m in obj_iter:
        tokens[m.group("object_name")] = m.group("value")
    """
    for k,v in tokens.items():
        print ('{} = {}'.format (k, v))
    """

    # Format the metadata tokens as JSON fields
    s = []
    s.append(MetadataSpecification)
    s.append(format_CollectionReference (tokens["SHORTNAME"], tokens["VERSIONID"]))
    s.append(format_DataGranule (tokens["LOCALGRANULEID"], tokens["PRODUCTIONDATETIME"]))
    s.append(format_GranuleUR (tokens["LOCALGRANULEID"]))
    s.append(format_InputGranules (tokens["INPUTPOINTER"]))
    s.append(format_PGEVersionClass (tokens["PGEVERSION"]))

    # Omit the SpatialExtent field when the input lacks a bounding polygon.
    if "GRINGPOINTLATITUDE" in tokens.keys():
        s.append(format_SpatialExtent (tokens["GRINGPOINTLATITUDE"], tokens["GRINGPOINTLONGITUDE"]))

    s.append(format_TemporalExtent (tokens["RANGEBEGINNINGDATE"], tokens["RANGEBEGINNINGTIME"],
                                    tokens["RANGEENDINGDATE"], tokens["RANGEENDINGTIME"]))

    # Join JSON fields into a single string
    s = "{%s}" % (",".join(s))

    # Validate the JSON formatted string
    parsed = json.loads(s)

    # Pretty-print the JSON output file
    s_pretty = json.dumps(parsed, indent=2)
    with open (jsonfile, "w") as fp:
        fp.write(s_pretty)

    return jsonfile

def filter_file_list (infile, outfile):

    if outfile is None:
        outfile = infile + ".out"

    # Read the input file list
    with open (infile, "r") as fp:
        infile_list = fp.read()
    infile_list = infile_list.split('\n')
    # filter out empty strings
    infile_list = list(filter (None, infile_list))

    # Write out the same file list, inserting a .cmr.json filename after each .met file
    with open (outfile, "w") as fp:
        for filename in infile_list:
            fp.write(filename + "\n")
            if filename.endswith (".met"):
                jsonfile = convert_odl_file (filename)
                fp.write (jsonfile + "\n")

def main():
    parser = argparse.ArgumentParser(description='Convert TEMPO ODL metadata to UMM-G JSON')
    parser.add_argument('--filter', default=None,
                        help="input filename list (--output is a filtered list, plus JSON files)")
    parser.add_argument('--input', default=None,
                        help="input ODL file (--output is one JSON file)")
    parser.add_argument('--output', default=None,
                        help="output file")
    args = parser.parse_args()

    if args.input is not None:
        convert_odl_file (args.input, jsonfile=args.output)
    elif args.filter is not None:
        filter_file_list (args.filter, args.output)
    else:
        parser.print_usage(sys.stderr)

if __name__ == "__main__":
    main()
