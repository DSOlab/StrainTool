#! /usr/bin/python2.7

############################################## standard libs
import sys
from copy import deepcopy
from math import degrees, radians, floor, ceil
##############################################  numpy & argparse
import numpy
import argparse
##############################################  pystrain
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *
############################################## ploting
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description='Test routine for z-weights.',
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
print '[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell))

##  Make a new station list (copy of the original one), where all coordinates
##+ are in UTM. All points should belong to the same ZONE.
mean_lon = degrees(sum([ x.lon for x in sta_list_ell ])/len(sta_list_ell))
utm_zone = floor(mean_lon/6)+31
utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
print '[DEBUG] Mean longtitude is {} deg.; using Zone = {} for UTM'.format(mean_lon, utm_zone)
sta_list_utm = deepcopy(sta_list_ell)
for idx, sta in enumerate(sta_list_utm):
    N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), utm_zone)
    sta_list_utm[idx].lon = E
    sta_list_utm[idx].lat = N
    assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
print '[DEBUG] Station list transformed to UTM.'

##  Get the barycenter of the points (in UTM and geodetic)
clon = degrees(sum([ x.lon for x in sta_list_ell ])/len(sta_list_ell))
clat = degrees(sum([ x.lat for x in sta_list_ell ])/len(sta_list_ell))
cN, cE, _, _ = ell2utm(radians(clat), radians(clon), Ellipsoid("wgs84"), utm_zone)
assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."

##  Call the z-weights function
zweights = z_weights(sta_list_utm, cE, cN, True)

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
my_map.plot(clon, clat, markersize=12)
for sta in sta_list_ell:
    x, y = my_map(degrees(sta.lon), degrees(sta.lat))
    my_map.plot(x, y, 'bo', markersize=10)
    plt.text(x, y, sta.name)
    my_map.drawgreatcircle(clon,clat,degrees(sta.lon),degrees(sta.lat), linewidth=2,color='g')
plt.show()
