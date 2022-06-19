#! /usr/bin/python3

import sys, os
import csv
import numpy as np
import datetime as dt

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as md

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

def read_safe_limits (filename):
    with open (filename) as csv_file:
        reader = csv.DictReader (csv_file)
        beg_sza = []
        beg_timet = []
        end_sza = []
        end_timet = []
        for row in list(reader)[1:]:
            beg_sza.append(float(row['beg_SZA']))
            beg_timet.append(float(row['beg_timet']))
            end_sza.append(float(row['end_SZA']))
            end_timet.append(float(row['end_timet']))

    beg = {'sza':np.asarray(beg_sza), 'timet':np.asarray(beg_timet)}
    end = {'sza':np.asarray(end_sza), 'timet':np.asarray(end_timet)}

    return beg, end

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Plot SZA at safe limit times')
    parser.add_argument('infile', help="CSV data file name")
    parser.add_argument('--label', help="plot title label")
    parser.add_argument('--outfile', help="plot file name")
    if len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    (beg, end) = read_safe_limits (args.infile)

    pdf = PdfPages (args.outfile)
    fig, ax1 = plt.subplots (1,1)
    fig = plt.figure(1)
    
    ax1.set_ylim (70.0, 160.0)
    ax1.set_xlabel ('Month')
    ax1.set_ylabel ('SZA [deg]')

    title = 'SZA @ TEMPO safety limits'
    if args.label is not None:
        title = '{} {}'.format(args.label,title)

    ax1.set_title (title)

    #xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')
    xfmt = md.DateFormatter('%b')
    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_major_locator(md.MonthLocator())

    #plt.subplots_adjust(bottom=0.2)
    #plt.xticks(rotation=25)

    beg_timet = beg["timet"]
    dates=[dt.datetime.fromtimestamp(ts) for ts in beg_timet]
    datenums=md.date2num(dates)
    ax1.plot(datenums, beg["sza"], label='morning')

    end_timet = end["timet"]
    dates=[dt.datetime.fromtimestamp(ts) for ts in end_timet]
    datenums=md.date2num(dates)
    ax1.plot(datenums, end["sza"], label='evening')
    
    ax1.legend(loc='lower left')
    
    pdf.savefig(fig)
    pdf.close()

if __name__ == '__main__':
    main()

