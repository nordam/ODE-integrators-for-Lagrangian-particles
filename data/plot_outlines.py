#!/usr/bin/env python

import sys
import netCDF4 as nc
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from datetime import datetime, timedelta
from matplotlib.patches import Polygon


def plot_boundary(datafilename, m, ax, color = None, label = None):
    with nc.Dataset(datafilename) as d:
        # Try some different common options for longitude
        # and latitude variable names.
        try:
            lons = d.variables['lon'][:,:]
            lats = d.variables['lat'][:,:]
        except KeyError:
            try:
                lons = d.variables['longitude'][:,:]
                lats = d.variables['latitude'][:,:]
            except KeyError:
                try:
                    lons = d.variables['gridLons'][:,:]
                    lats = d.variables['gridLats'][:,:]
                except KeyError:
                    lons = d.variables['lon_u'][:,:]
                    lats = d.variables['lat_u'][:,:]

        # Create list of points to draw outlines
        x, y = [], []
        x += list(lons[0,:])
        y += list(lats[0,:])
        x += list(lons[:,-1])
        y += list(lats[:,-1])
        x += list(lons[-1,:][::-1])
        y += list(lats[-1,:][::-1])
        x += list(lons[:,0][::-1])
        y += list(lats[:,0][::-1])
        X, Y = m(np.array(x), np.array(y))
        # Plot the outline
        if color is not None:
            m.plot(X, Y, label = label, color = color)
        else:
            m.plot(X, Y, label = label)


def main():
    # setup stereographic basemap.
    # lat_ts is latitude of true scale.
    # lon_0,lat_0 is central point.
    fig = plt.figure(figsize = (6,4))
    ax  = fig.add_subplot(111)
    width = 1500000
    m = Basemap(width=width,height=0.7*width,
            resolution='i',projection='stere',
            lat_ts=90,lat_0=60,lon_0=2, ax = ax)

    seacolor  = '#47C9FF'
    landcolor = '#9DDDC7'
    linewidth = 0.1
    alpha     = 0.3

    m.drawmapboundary(fill_color=seacolor, linewidth = 1)
    m.fillcontinents(color=landcolor, lake_color=seacolor, alpha = 1)
    m.drawcoastlines(linewidth = linewidth)
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,1.), linewidth = linewidth, alpha = alpha, labels = [1,0,0,0])
    m.drawmeridians(np.arange(-180.,181.,2.), latmax = 90, linewidth = linewidth, alpha = alpha, labels = [0,0,0,1])

    files = OrderedDict({})
    files['Selection from NorKyst 800 m'] = 'norkyst800.nc'
    files['Selection from Nordic 4 km'] = 'nordic4.nc'
    files['Selection from Arctic 20 km'] = 'arctic20.nc'
    for k, v in files.items():
        plot_boundary(v, m, ax, label = k)

    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig('outlines.png', dpi = 180)


if __name__ == '__main__':
    main()
