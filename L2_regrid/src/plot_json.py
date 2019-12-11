#! /usr/bin/env python

import sys
import os

from mpl_toolkits.basemap import Basemap

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

import numpy as np
import json

import argparse
parser = argparse.ArgumentParser(description='plot polygon vertices')
parser.add_argument('--format', default='pdf',
                                help="plot format [e.g. pdf, svgz]")
parser.add_argument('--lonlat', default=False, action='store_true',
                                help="plot in longitude-latitude coordinates")
parser.add_argument('filename')
args = parser.parse_args()

filename = args.filename
basename = os.path.basename(filename)

with open (filename) as f:
    j = json.load(f)

geom = j['geometries']
idlist = [geom[k]['id'] for k in range(len(geom))]

fig, ax = plt.subplots()
ax.ticklabel_format (style='sci', scilimits=(-3,4), useMathText=True)
ax.set_title(basename)  # <-- error message 'Unable to parse the pattern' is a known python bug

if args.lonlat:
    # plot lon-lat coordinates in degrees
    ax.set_xlabel ('longitude [deg]')
    ax.set_ylabel ('latitude [deg]')

    meters_per_mile = 5280.0 * 12.0 / 2.54 / 100.0
    width_miles = 100.0
    height_miles = 100.0
    m = Basemap(width=width_miles*meters_per_mile,
        height=height_miles*meters_per_mile,
        resolution='h',projection='aea',
        lat_1=29.5, lat_2=45.5, lat_0=37.5, lon_0=-96.0)
else:
    # plot Albers coordinates in km
    ax.set_xlabel ('x [km]')
    ax.set_ylabel ('y [km]')

patches = []
for g in geom[:]:
    coords = np.array(g['coordinates'][0])

    if args.lonlat:
        # plot lon-lat coordinates in degrees
        lon,lat = m(coords[:,0], coords[:,1], inverse=True)
        p_coords = zip(lon,lat)
    else:
        # plot Albers coordinates in km
        p_coords = np.array(coords) / 1.e3  # [km]

    poly = Polygon(p_coords, True, edgecolor='k')
    patches.append (poly)

p = PatchCollection (patches, alpha=0.4, linewidth=0.125)
p.set_array (np.array(idlist))

ax.add_collection(p)
ax.autoscale_view()

#plt.show()
fig.savefig (basename + "." + args.format, bbox_inches='tight')

