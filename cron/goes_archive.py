#! /usr/bin/env python3

# for eprint definition
from __future__ import print_function

import os, sys
import glob
import hashlib
import re
import pathlib
from datetime import datetime,timedelta
import argparse
import subprocess

from netCDF4 import Dataset as NetCDFFile

Verbose = False
Dryrun = False

# python3 will provide file= redirection to stderr
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def suffix_sha1 (path):
    return path + ".sha1"

def sha1(fname):
    hash_sha1 = hashlib.sha1()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha1.update(chunk)
    return hash_sha1.hexdigest()

def have_valid_checksum (ncfile_path):
    if not os.path.isfile(ncfile_path):
        eprint ("*** Error: invalid filename path: {}".format(ncfile_path))
        return False
    ncfile_sha1sum = sha1 (ncfile_path)

    ncfile_dir = os.path.dirname (ncfile_path)
    ncfile_basename = os.path.basename (ncfile_path)
    chksum_path = os.path.join (ncfile_dir, suffix_sha1(ncfile_basename))

    if os.path.isfile (chksum_path):
        with open (chksum_path, "r") as f:
            read_sha1 = f.readline().strip()
        result = (read_sha1 == ncfile_sha1sum)
        if result is False:
            eprint ("*** bad checksum: {} {}".format(ncfile_sha1sum, ncfile_path))
        return result

    try:
        print ("creating {}".format(chksum_path))
        if not Dryrun:
            with open (chksum_path, "w") as f:
                f.write (ncfile_sha1sum)
            return True
    except:
        eprint ("*** Error: creating {}".format(chksum_path))
        return False

def move_file_to_dir (path, dest_dir, clobber=True):
    new_path = os.path.join (dest_dir, os.path.basename(path))
    if os.path.isfile(new_path) and not clobber:
        eprint ("file exists: {}".format(new_path))
        return
    if not Dryrun:
        os.rename (path, new_path)
    if Verbose:
        print ("{} <- {}".format(dest_dir, path))

def get_yday_subdir (base, zone):
    """
    It's convenient to have one directory contain all GOES imagery
    needed for a single operational day.  To implement this, we organize
    the files using the day-of-year in the TEMPO satellite-local time zone.
    """
    # parse filename to get timestamp to whole seconds precision
    filename_regex = 'OR_ABI-L2-CMIPF-M\d{1}C\d{2}_G\d{2}_s(\d{13})\d_'
    fields = re.search (filename_regex, base)
    if fields is None:
        eprint ("*** Error: regex mismatch: {}".format(base))
        return None
    # help strptime by adding a timezone specification
    tstamp = fields.group(1)+'+0000'

    # archive using day-of-year in satellite local time zone (SLT)
    tstamp_obj = datetime.strptime(tstamp, '%Y%j%H%M%S%z')
    tstamp_slt_obj = tstamp_obj - timedelta(hours=zone)
    yday_subdir = datetime.strftime (tstamp_slt_obj, "%Y/%j")
    return yday_subdir

def classify_goes_file (path, root_dir):
    try:
        with NetCDFFile(path, "r") as nc:
            orbital_slot = nc.getncattr ("orbital_slot")
    except:
        eprint ("*** Error: cannot read orbital slot from file: {}".format(path))
        return None, None

    if orbital_slot == "GOES-East":
        goes_subdir = "east_cmi"
        dbfilename  = "cmieast.sqlite"
    elif orbital_slot == "GOES-West":
        goes_subdir = "west_cmi"
        dbfilename  = "cmiwest.sqlite"
    else:
        eprint ("*** Error: unsupported orbital slot {}:{}".format(orbital_slot, path))
        return None, None

    return goes_subdir, os.path.join (root_dir, dbfilename)

def archive_goes_file (path, root_dir, zone):
    """
    Classify the file, decide where it goes, move it there,
    and register the final archived path in the tracking database
    """
    # classify the file as either GOES-East or GOES-West
    (goes_subdir, dbfile) = classify_goes_file (path, root_dir)
    if goes_subdir is None:
        return

    # Determine the archive path from the filename
    base = os.path.basename (path)
    yday_subdir = get_yday_subdir (base, zone)
    if yday_subdir is None:
        return

    # Ensure destination directory exists
    dest_dir = os.path.join (root_dir, yday_subdir, goes_subdir)
    if not Dryrun:
        pathlib.Path(dest_dir).mkdir(parents=True, exist_ok=True)

    # Move files into the destination directory
    move_file_to_dir (path, dest_dir, clobber=False)
    move_file_to_dir (suffix_sha1(path), dest_dir, clobber=False)

    # register final file path in sqlite database
    final_path = os.path.realpath (os.path.join (dest_dir, base))
    argv = ["asdc_files.py", "--dbfile", dbfile, "--add", final_path]
    if not Dryrun:
        result = subprocess.run (argv, check=True)
    else:
        print (" ".join (argv))

def handle_bad_checksum (path, root_dir):
    dir = os.path.join (root_dir, "bad_checksum")
    pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
    move_file_to_dir (path, dir, clobber=False)
    move_file_to_dir (suffix_sha1(path), dir, clobber=False)

def main():
    default_zone = 6 # spacecraft local time zone
    parser = argparse.ArgumentParser(description='Validate and archive GOES CMI data products')
    parser.add_argument('--root', metavar='DIR', help="Root directory to hold GOES archive directory tree: YYYY/ddd/xxxx_cmi")
    parser.add_argument('--zone', type=int, default=default_zone, help="Satellite local time zone [integer hours east of UTC; default={}]".format(default_zone))
    parser.add_argument('--verbose', help="Show verbose progress", action='store_true')
    parser.add_argument('--dryrun', help="Print actions, but don't modify any files or directories", action='store_true')
    parser.add_argument('source', metavar='SOURCE', help="source directory or source file list", nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    root_dir = args.root
    if not os.path.isdir (root_dir):
        eprint ("*** Error: cannot access destination root directory: {}".format(root_dir))
        sys.exit(1)

    global Verbose
    Verbose = args.verbose

    global Dryrun
    Dryrun = args.dryrun

    if Dryrun:
        Verbose = True

    if len(args.source) == 1 and os.path.isdir(args.source[0]):
        ncfile_list = glob.glob (os.path.join (args.source[0], "*.nc"))
    else:
        ncfile_list = args.source

    num_bad_checksums = 0
    for path in ncfile_list:
        if have_valid_checksum (path):
            archive_goes_file (path, root_dir, args.zone)
        else:
            handle_bad_checksum (path, root_dir)
            num_bad_checksums += 1

    if num_bad_checksums > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()
