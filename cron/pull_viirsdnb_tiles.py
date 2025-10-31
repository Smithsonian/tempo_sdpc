#!/usr/bin/env python3

# Attempts to do HTTP Gets with urllib2(py2) urllib.requets(py3) or subprocess
# if tlsv1.1+ isn't supported by the python ssl module
#
# Will download csv or json depending on which python module is available
#

#from __future__ import (division, print_function, absolute_import, unicode_literals)

import argparse
import os
import os.path
import shutil
import sys
import time
import re

from io import StringIO

DESC = "This script will download selected files and store them to the specified path"

USERAGENT = 'Custom user agent'

# this is the choice of last resort, when other attempts have failed
def getcURL(url, headers=None, out=None):
    # OS X Python 2 and 3 don't support tlsv1.1+ therefore... cURL
    import subprocess
    try:
        print('trying cURL', file=sys.stderr)
        args = ['curl', '--fail', '--retry', '5', '--retry-connrefused', '-sS', '-L', '-b session', '--get', url]
        for (k,v) in headers.items():
            args.extend(['-H', ': '.join([k, v])])
        if out is None:
            # python3's subprocess.check_output returns stdout as a byte string
            result = subprocess.check_output(args)
            return result.decode('utf-8') if isinstance(result, bytes) else result
        else:
            subprocess.call(args, stdout=out)
    except subprocess.CalledProcessError as e:
        print('curl GET error message: %' + (e.message if hasattr(e, 'message') else e.output), file=sys.stderr)
    return None

# read the specified URL and output to a file
def geturl(url, token=None, out=None):
    headers = { 'user-agent' : USERAGENT }
    if not token is None:
        headers['Authorization'] = 'Bearer ' + token
    try:
        import ssl
        CTX = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
        if sys.version_info.major == 2:
            import urllib2
            try:
                fh = urllib2.urlopen(urllib2.Request(url, headers=headers), context=CTX)
                if out is None:
                    return fh.read()
                else:
                    shutil.copyfileobj(fh, out)
            except urllib2.HTTPError as e:
                print('TLSv1_2 sys 2 : HTTP GET error code: %d' % e.code, file=sys.stderr)
                return getcURL(url, headers, out)
            except urllib2.URLError as e:
                print('TLSv1_2 sys 2 : Failed to make request: %s, RETRYING' % e.reason, file=sys.stderr)
                return getcURL(url, headers, out)
            return None

        else:
            from urllib.request import urlopen, Request, URLError, HTTPError
            try:
                fh = urlopen(Request(url, headers=headers), context=CTX)
                if out is None:
                    return fh.read().decode('utf-8')
                else:
                    shutil.copyfileobj(fh, out)
            except HTTPError as e:
                print('TLSv1_2 : HTTP GET error code: %d' % e.code, file=sys.stderr)
                return getcURL(url, headers, out)
            except URLError as e:
                print('TLSv1_2 : Failed to make request: %s' % e.reason, file=sys.stderr)
                return getcURL(url, headers, out)
            return None

    except AttributeError:
      return getcURL(url, headers, out)

Parse_Name = re.compile (r'VNP46A3\.A\d{7}\.h(?P<h>\d{2})v(?P<v>\d{2}).002\.\d{13}\.h5')

def is_selected_file (fname):
    m = Parse_Name.match (fname)
    if m is None:
        return False
    h = int(m.group('h'))
    v = int(m.group('v'))
    return 4 <= h and h <= 13 and 2 <= v and v <= 7

def sync(src, dest, dryrun, tok):
    '''synchronize src url with dest directory'''
    try:
        import csv
        files = {}
        files['content'] = [ f for f in csv.DictReader(StringIO(geturl('%s.csv' % src, tok)), skipinitialspace=True) ]
    except ImportError:
        import json
        files = json.loads(geturl(src + '.json', tok))

    num_download_files = 0
    total_download_size_MB = 0

    for f in files['content']:
        fname = f['name']
        if not is_selected_file (fname):
            continue
        path = os.path.join(dest, fname)
        filesize = int(f['size'])
        total_download_size_MB += filesize / 1.e6
        num_download_files += 1
        if dryrun:
            print ('%7d,%s' % (filesize, path))
            continue
        url = src + '/' + fname
        if filesize == 0:
            # filesize=0 from response indicates a directory
            print('skipping zero size file: %s' % (fname))
        else:
            try:
                if not os.path.exists(path) or os.path.getsize(path) == 0:
                    print('downloading: ' , path)
                    with open(path, 'w+b') as fh:
                        geturl(url, tok, fh)
                else:
                    print('skipping: ', path)
            except IOError as e:
                print("open `%s': %s" % (e.filename, e.strerror), file=sys.stderr)
                sys.exit(-1)

    if dryrun:
        print('Would download %d files, total size = %ld MB' % (num_download_files, total_download_size_MB))
    return 0

def read_token (token_file):
    if not os.path.isfile (token_file):
        print('nonexistent file: %s' % (token_file))
        return None
    with open (token_file, "r") as fp:
        lines = fp.readlines()
    lines = [s.strip('\n') for s in lines]
    return "".join (lines)

def _main(argv):
    """
    For example, source:
        https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/VNP46A3/2025/213/
    and token:
        $HOME/.ssh/LAADS_EDL_token.txt
    """
    parser = argparse.ArgumentParser(prog=argv[0], description=DESC)
    parser.add_argument('-s', '--source', dest='source', metavar='URL', help='Recursively download files at URL', required=True)
    parser.add_argument('-d', '--destination', dest='destination', metavar='DIR', help='Store directory structure in DIR', required=True)
    parser.add_argument('-t', '--token', dest='token', metavar='TOKFILE', help='Use app token file TOKFILE to authenticate', required=True)
    parser.add_argument('--dryrun', action='store_true', help='Only print file size and destination, do not download')
    args = parser.parse_args(argv[1:])

    if not os.path.exists(args.destination):
        os.makedirs(args.destination)

    the_token = read_token (args.token)
    if the_token is None:
        return -1

    return sync(args.source, args.destination, args.dryrun, the_token)

if __name__ == '__main__':
    try:
        sys.exit(_main(sys.argv))
    except KeyboardInterrupt:
        sys.exit(-1)
