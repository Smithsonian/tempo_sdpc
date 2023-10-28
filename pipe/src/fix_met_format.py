#! /usr/bin/env python3

import sys
import re
import argparse

def impose_paren_state (name, parens_required, text):
    """
    Impose a specified paren state on the ODL value for variable 'name'.
    Some ODL values are defined as scalar and will fail if parens are used.
    Other ODL values are defined as array and will fail if parens are not used.
    The two possible states are parens_required and parens_illegal.
    If parens are missing, but needed, insert them.
    If parens are present, but not wanted, remove them.
    Otherwise, do nothing.
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
    # If a '$' is present, the field has an unexpanded macro, so nothing
    if s.find('$') >= 0:
        return text
    if parens_required:
        s = s.strip()
        r = text[b:eol]
        if s[0] != '(':
            r = ' (' + r
        if s[-1] != ')':
            r = r + ' )'
        s = r
    else:
        if s.find(')') >= 0:
            s = s.rsplit(')',1)[0]
        if s.find('(') >= 0:
            s = " " + s.split('(',1)[1]

    return text[:b] + s + text[eol:]

def filter_boundary_coords_for_asdc (name, text):
    """
    ASDC asked us to remove newlines from bounding polygon coordinate arrays.
    Whatever.
    """
    # find the object
    beg = re.search ('OBJECT[ ]*=[ ]*{}'.format(name), text)
    if (beg == None):
        return text
    b = beg.end()
    end = re.search ('[ \n]*END_OBJECT[ ]*=[ ]*{}'.format(name), text)
    e = end.start()
    # find the value
    value = re.search ('VALUE[ ]*=', text[b:e])
    b += value.end()
    s = " " + " ".join(text[b:e].split())
    return text[:b] + s + text[e:]

def process_file (metfile):

    with open (metfile, 'r') as fp:
        text = fp.read()

    # True means parens required -- insert parens if they don't exist
    # False means parens illegal -- remove parens if they exist
    paren_states = {'VERSIONID': False, 'PGEVERSION': False, 'INPUTPOINTER': True}

    for name,state in paren_states.items():
        text = impose_paren_state (name, state, text)

    coordinates = ['GRINGPOINTLATITUDE', 'GRINGPOINTLONGITUDE', 'GRINGPOINTSEQUENCENO']
    for name in coordinates:
        text = filter_boundary_coords_for_asdc (name, text)

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
