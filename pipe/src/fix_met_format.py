#! /usr/bin/env python3

import sys
import re
import argparse

def fix_parens (name, text):
    """
    Remove parentheses from ODL scalar value for variable 'name'
    """
    # find the object
    beg = re.search ('OBJECT[ ]*=[ ]*{}'.format(name), text)
    if (beg == None):
        return text
    end = re.search ('END_OBJECT[ ]*=[ ]*{}'.format(name), text)
    b = beg.end()
    e = end.start()
    # find the value
    value = re.search ('VALUE[ ]*=', text[b:e])
    b += value.end()
    eol = b + text[b:e].find('\n')
    s = text[b:eol]
    if s.find(')') >= 0:
        s = s.rsplit(')',1)[0]
    if s.find('(') >= 0:
        s = " " + s.split('(',1)[1]
    return text[:b] + s + text[eol:]

def process_file (metfile):

    with open (metfile, 'r') as fp:
        text = fp.read()

    name_list = {'VERSIONID', 'PGEVERSION'}

    for n in name_list:
        text = fix_parens (n, text)

    with open (metfile, 'w') as fp:
        fp.write(text)

def main ():
    parser = argparse.ArgumentParser(description='fix SDPTK .met formatting bugs')
    parser.add_argument('metfile', help="MET file name")
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()
    process_file (args.metfile)

if __name__ == '__main__':
    main()
