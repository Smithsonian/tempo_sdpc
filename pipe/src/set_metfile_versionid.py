#! /usr/bin/env python3

import sys
import re
import argparse

def process_file (metfile, versionid):

    with open (metfile, 'r') as fp:
        text = fp.read()

    beg = re.search (r'=\s+VERSIONID', text)
    b = beg.start()
    d = re.match (r'=\s+VERSIONID\s+NUM_VAL\s+=\s+1\s+VALUE\s+=\s+(\d+)', text[b:])

    with open (metfile, 'w') as fp:
        fp.write(text[:b+d.start(1)] + versionid + text[b+d.end(1):])

def main():
    parser = argparse.ArgumentParser(description='Set .met VERSIONID entry')
    parser.add_argument('--versionid', help="VERSIONID value")
    parser.add_argument('metfile', help="MET file name")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()
    process_file (args.metfile, args.versionid)

if __name__ == '__main__':
    main()
