#! /usr/bin/python2.7

from __future__ import print_function
import sys
from copy import deepcopy
from math import degrees, radians, floor, ceil
import numpy
import argparse
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *
from mpl_toolkits.basemap import Basemap

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description='Test routine for delauney triangulation.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
Send bug reports to:
  Xanthos Papanikolaou, xanthos@mail.ntua.gr
  Demitris Anastasiou,danast@mail.ntua.gr
November, 2017'''))

parser.add_argument('-i', '--input-file',
    default=None,
    metavar='INPUT_FILE',
    dest='gps_file',
    required=True)

##  Parse command line arguments.
args = parser.parse_args()

##  Parse stations from input file
sta_list_ell = parse_ascii_input(args.gps_file)
print('[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell)))

##  Make a new station list (copy of the original one), where all coordinates
##+ are in UTM. All points should belong to the same ZONE.
mean_lon = degrees(sum([ x.lon for x in sta_list_ell ])/len(sta_list_ell))
utm_zone = floor(mean_lon/6)+31
utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
print('[DEBUG] Mean longtitude is {:+8.4f} deg.; using Zone = {:+8.4f} for UTM'.format(mean_lon, utm_zone))
sta_list_utm = deepcopy(sta_list_ell)
for idx, sta in enumerate(sta_list_utm):
    N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), utm_zone)
    sta_list_utm[idx].lon = E
    sta_list_utm[idx].lat = N
    assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
print('[DEBUG] Station list transformed to UTM.')

##  Get the barycenter of the points (in UTM and geodetic)
clon = degrees(sum([ x.lon for x in sta_list_ell ])/len(sta_list_ell))
clat = degrees(sum([ x.lat for x in sta_list_ell ])/len(sta_list_ell))
cN, cE, _, _ = ell2utm(radians(clat), radians(clon), Ellipsoid("wgs84"), utm_zone)
assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."

points = numpy.array([ [sta.lon, sta.lat] for sta in sta_list_utm ])
tri = Delaunay(points)
print('Here are the simplices (i.e. indices of the points forming the simplices in the triangulation.')
for idx, trng in enumerate(tri.simplices):
    print('idx:{:2d} triangle indexes {:2d}->{:2d}->{:2d} (aka {:}-{:}-{:})'.format(idx, trng[0], trng[1], trng[2], sta_list_ell[trng[0]].name, sta_list_ell[trng[1]].name, sta_list_ell[trng[2]].name))
print('\nHere are the coordinates of each triangle:')
for idx, trng in enumerate(points[tri.simplices]):
    print('idx:{:2d} triangle coordinates: ({:10.6f}, {:10.6f})->({:10.6f}, {:10.6f})->({:10.6f}, {:10.6f})'.format(idx, trng[0][0], trng[0][1], trng[1][0], trng[1][1], trng[2][0], trng[2][1],))
print('\nHere are the coordinates of each point (validation):')
for s in sta_list_utm:
    print('{:} {:10.6f} {:10.6f}'.format(s.name, s.lon, s.lat))

##  Plot
lat0    = degrees(sum([ x.lat for x in sta_list_ell ])/len(sta_list_ell))
lon0    = degrees(sum([ x.lon for x in sta_list_ell ])/len(sta_list_ell))
lons    = [ degrees(x.lon) for x in sta_list_ell ]
lats    = [ degrees(x.lat) for x in sta_list_ell ]
lon_off = (max(lons)-min(lons))/10
lat_off = (max(lats)-min(lats))/10
my_map = Basemap(projection='merc', lat_0 = lat0, lon_0 = lon0, resolution = 'c', llcrnrlon=min(lons)-lon_off, llcrnrlat=min(lats)-lat_off, urcrnrlon=max(lons)+lon_off, urcrnrlat=max(lats)+lat_off)
my_map.drawcoastlines()
my_map.drawcountries()
my_map.fillcontinents(color = 'coral')
my_map.drawmapboundary()
my_map.drawmeridians(numpy.arange(floor(min(lons)), ceil(max(lons)), 2), labels=[True,False,False,True])
my_map.drawparallels(numpy.arange(floor(min(lats)), ceil(max(lats)), 2), labels=[False,True,True,False], fontsize=10)
for sta in sta_list_ell:
    x, y = my_map(degrees(sta.lon), degrees(sta.lat))
    my_map.plot(x, y, 'bo', markersize=10)
    plt.text(x, y, sta.name)
pts_ell_map = numpy.array([ list(my_map(degrees(sta.lon), degrees(sta.lat))) for sta in sta_list_ell ])
plt.triplot(pts_ell_map[:,0], pts_ell_map[:,1], tri.simplices.copy())
plt.show()
