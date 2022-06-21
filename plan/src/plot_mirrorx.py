#! /usr/bin/python3

import sys, os
import csv
import numpy as np
import datetime as dt
from dateutil import tz, parser
import re

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as md

matplotlib.rcParams['hatch.linewidth'] = 0.2

# Axis label numbers are still serif font.  Why?
plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{lmodern}\renewcommand*\familydefault{\sfdefault}\usepackage[T1]{fontenc}')

def read_plan (filename):
    with open (filename) as csv_file:
        reader = csv.DictReader (filter(lambda row: row[0]!='#', csv_file))
        label = []
        time = []
        duration = []
        mirror_x = []
        num_steps = []
        repeat = []
        timestamp = []
        for row in list(reader):
            label.append(int(row['label']))
            time.append(float(row['time']))
            duration.append(float(row['duration']))
            mirror_x.append(float(row['mirror_x']))
            num_steps.append(int(row['num_steps']))
            repeat.append(int(row['repeat']))
            timestamp.append(row['timestamp'])

    plan_id = None
    with open (filename, 'r') as fp:
        for line in fp:
            m = re.search(r"plan_id\W=\W(\w+)", line)
            if m is not None:
                plan_id = m.group(1)
                break

    plan = {'label':np.asarray(label),
            'time':np.asarray(time),
            'duration':np.asarray(duration),
            'mirror_x':np.asarray(mirror_x),
            'num_steps':np.asarray(num_steps),
            'repeat':np.asarray(repeat),
            'timestamp':np.asarray(timestamp),
            'plan_id':plan_id}

    return plan

def plot_scan(ax, plan, step_size, indices):
    xs = plan['mirror_x'][indices]/1000.0
    t0s = plan['time'][indices]
    dts = plan['duration'][indices]
    ns = plan['num_steps'][indices]
    nrs = plan['repeat'][indices]

    t0s = t0s - t0s[0]
    xf = xs - step_size * ns

    ax.set_xlabel (r'$\delta t$ [hour]')
    ax.set_yticks ([-25, 0, 25])
    ax.text (0.01, 0.04, plan['timestamp'][indices[0]], fontsize='xx-small', transform=ax.transAxes)

    # plot scans
    for k in indices - indices[0]:
        t0 = t0s[k]/3600.0
        dt = dts[k]/3600.0
        x0 = xs[k]
        x1 = xf[k]
        if nrs[k] == 0:
            ax.plot ([t0, t0+dt], [x0, x1], color='k', linewidth=0.5)
        else:
            for i in np.arange(nrs[k]):
                ax.plot ([t0, t0+dt], [x0, x1], color='b', linewidth=0.5)
                t0 += dt

    ymin, ymax = ax.get_ylim()

    # highlight gaps
    t1s = t0s + dts * np.where(nrs > 0, nrs, 1)
    gaps = np.argwhere (t0s[1:] - t1s[:-1])

    for j in gaps:
        t0 = t1s[j][0]/3600.0
        t1 = t0s[j+1][0]/3600.0
        ax.fill_between ([t0, t1], [ymin, ymin], y2=[ymax,ymax], hatch='////',
                         linewidth=0.2, color='r', fc='w')

def new_page (plan_id):
    fig, axs = plt.subplots (7,1, sharex=True)
    fig.subplots_adjust (hspace=0)
    fig.text(0.04, 0.5, r'\verb|mirror_x| [mrad]', va='center', rotation='vertical')
    fig.text(0.125, 0.9, r'\verb|plan_id|={}'.format(plan_id), va='center')
    return fig, axs

def find_yday (timestamps, zone_name):
    to_zone = tz.gettz(zone_name)

    doy = []
    for ts in timestamps:
        utc = parser.isoparse(ts)
        t_zone = utc.astimezone (to_zone)
        tt = t_zone.timetuple()
        doy.append(tt.tm_yday)

    return doy

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Plot scan mirror_x vs time')
    parser.add_argument('infile', help="CSV data file name")
    parser.add_argument('--step', type=float, default=0.057, help="step size [mrad]")
    parser.add_argument('--label', help="plot title label")
    parser.add_argument('--outfile', help="plot file name", default="mirrorx.pdf")
    if len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    plan = read_plan (args.infile)
    plan_id = plan['plan_id']

    ydays = find_yday (plan['timestamp'], 'America/Chicago')

    pdf = PdfPages (args.outfile)
    fig, axs = new_page(plan_id)

    k = 0
    for day in np.unique(ydays):
        indices = np.flatnonzero (ydays == day)
        plot_scan (axs[k], plan, args.step, indices)
        k = (k+1) % 7
        if k == 0:
            pdf.savefig(fig)
            fig, axs = new_page(plan_id)

    if k != 0:
        pdf.savefig(fig)
    pdf.close()

if __name__ == '__main__':
    main()

