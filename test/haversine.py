#! /usr/bin/python

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

print('Mean long/lat = {:8.3f}, {:8.3f}'.format(clon, clat))
brc = Station(lat=cN, lon=cE)
bre = Station(lon=radians(clon), lat=radians(clat))
dr_cart = [ x.distance_from(brc)[2] for x in sta_list_utm ]
dr_ellp = [ x.haversine_distance(bre)[0] for x in sta_list_ell ]

for idx, dr in enumerate(dr_cart):
    print('lon={:+10.5f} lat={:10.5f}'.format(degrees(sta_list_ell[idx].lon), degrees(sta_list_ell[idx].lat)))
    print('Cartesian distance {:10.1f}, on sphere: {:10.1f}, diff={:5.1f}km (negative lon {:})'.format(dr/1e3, dr_ellp[idx]/1e3, abs(dr-dr_ellp[idx])/1e3, sta_list_ell[idx].lon<0e0))
