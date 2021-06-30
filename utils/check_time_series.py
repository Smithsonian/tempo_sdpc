#! /usr/bin/env python

import os, sys
import glob
import argparse
import netCDF4
import numpy as np

def get_file_times (fname, varpath):
    dset = netCDF4.Dataset(fname, 'r')
    tstart = dset.time_coverage_start_since_epoch
    gpath = os.path.dirname (varpath)
    vname = os.path.basename (varpath)
    if len(gpath) > 1:
        g = dset[gpath]
    else:
        g = dset
    t = g.variables[vname][:]
    dset.close ()
    return tstart,t

def check_for_gaps (fprev, f, varpath, before):

    have_prev = not (fprev == None)

    tstart,t = get_file_times (f, varpath)
    if len(t) == 0:
        return []

    dt = np.diff(t)
    mean_dt = np.mean(dt)
    dt_thresh = 3 * mean_dt
    num_bad_dt = len(np.where(dt < 0)[0])
    where_big_dt = np.where(dt > dt_thresh)
    num_big_dt = len(where_big_dt[0])

    vb = os.path.basename (varpath)
    fb = os.path.basename (f)

    issues = []

    if have_prev:
        fbprev = os.path.basename (fprev)
        tstart_prev,tprev = get_file_times (fprev, varpath)
        if tstart < tstart_prev:
            issues.append ('file out of order')
        if len(tprev) > 1 and ((t[0] - tprev[-1]) > dt_thresh):
            issues.append ('%0.3f sec gap preceeding' % (t[0]-tprev[-1]))

    if t[0] + before > tstart:
        issues.append ('insufficient padding: tstart-t[0] = %0.3f' % (tstart-t[0]))
    if num_bad_dt > 0:
        issues.append ('%d non-monotonic entries' % (num_bad_dt))
    if num_big_dt > 0:
        issues.append ('%d gaps > %0.3f sec (max = %0.3f)' % (num_big_dt, dt_thresh, np.max(dt)))
        kv = where_big_dt[0]
        for k in kv:
            issues.append ('gap: %s[%d] = %0.3f' % (vb, k, t[k]))
            issues.append ('     %s[%d] = %0.3f' % (vb, k+1, t[k+1]))

    return issues

def main():
    parser = argparse.ArgumentParser(description='check time series for continuity')
    parser.add_argument('--varpath', help="path to time series variable")
    parser.add_argument('--before', type=float, default=0.0, help="expected padding [sec]")
    parser.add_argument('files', nargs=argparse.REMAINDER)
    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    n = len(args.files)
    if n < 1:
        parser.print_usage(sys.stderr)
        sys.exit(0)

    print ('Variable: {}'.format(args.varpath))

    num_ok = 0
    f_prev = None

    for f in args.files:
        issues = check_for_gaps (f_prev, f, args.varpath, args.before)
        if len(issues) > 0:
            if num_ok > 0:
                print ('{} previous files OK'.format(num_ok))
            num_ok = 0
            print (os.path.basename(f))
            print ('\n'.join(issues))
        else:
            num_ok = num_ok + 1
        f_prev = f

    if num_ok == n:
        print('OK')

if __name__ == "__main__":
    main()

