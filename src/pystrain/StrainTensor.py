#! /usr/bin/python2.7
#-*- coding: utf-8 -*-

from __future__ import print_function
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
import pystrain.grid
############################################## ploting
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

def cut_rectangle(xmin, xmax, ymin, ymax, sta_lst):
    new_sta_lst = []
    for sta in sta_lst:
        if sta.lon >= xmin and sta.lon <= xmax and sta.lat >= ymin and sta.lat <= ymax:
            new_sta_lst.append(sta)
    return new_sta_lst

def gmt_script(input_file, sta_lst, tensor_lst, utm_zone, outfile='gmt_script', projscale=9000000, strsc=50, frame=2):
    lons    = [ degrees(x.lon) for x in sta_lst ]
    lats    = [ degrees(x.lat) for x in sta_lst ]
    west    = min(lons)
    east    = max(lons)
    south   = min(lats)
    north   = max(lats)
    gmt_range=('-R{:5.2f}/{:5.2f}/{:5.2f}/{:5.2f}'.format(west, east, south, north)) 
    gmt_scale=('-Lf20/33.5/36:24/100+l+jr')
    gmt_proj=('-Jm24/37/1:{:}'.format(projscale))
    with open(outfile, 'w') as fout:
        print('#!/bin/bash', file=fout)
        print('outfile=strain_rates.ps', file=fout)
        print('GMTVRB=2', file=fout)
        print('VSC=\"0.05\"', file=fout)
        print('VHORIZONTAL=0', file=fout)
        print('LABELS=0', file=fout)
        print('while [ $# -gt 0 ] ; do case \"$1\" in', file=fout)
        print('-a)\n  VHORIZONTAL=1\n  shift\n  ;;', file=fout)
        print('-l)\n  LABELS=1\n  shift\n  ;;', file=fout)
        print('*)\n  shift\n  ;;', file=fout)
        print('esac\ndone', file=fout)
        print('gmt gmtset MAP_FRAME_TYPE fancy', file=fout)
        print('gmt gmtset PS_PAGE_ORIENTATION portrait', file=fout)
        print('gmt gmtset FONT_ANNOT_PRIMARY 10', file=fout)
        print('gmt gmtset FONT_LABEL 10', file=fout)
        print('gmt gmtset MAP_FRAME_WIDTH 0.12c', file=fout)
        print('gmt gmtset FONT_TITLE 18p,Palatino-BoldItalic', file=fout)
        print('gmt gmtset PS_MEDIA 22cx22c', file=fout)
        print('gmt psbasemap {:} {:} {:} -B{:}:."Tensors": -P -K > $outfile'.format(gmt_range, gmt_proj, gmt_scale, frame), file=fout)
        print('gmt pscoast -R -J -O -K -W0.25 -G195 -Df -Na >> $outfile', file=fout)
        print('awk \'{{print $3,$2}}\' .station.info.dat | gmt psxy -Jm -O -R -Sc0.10c -W0.005c -Ggreen -K >>$outfile', file=fout)
        print("awk \'{{print $3,$2,0,$6,$8+90}}\' .strain.info.dat | gmt psvelo -Jm {:} -Sx{:} -L -A10p+e -Gblue -W2p,blue -V${{GMTVRB}} -K -O>> $outfile".format(gmt_range, strsc), file=fout)
        print('awk \'{{print $3,$2,$4,0,$8+90}}\' .strain.info.dat | gmt psvelo -Jm {:} -Sx{:} -L -A10p+e -Gred -W2p,red -V${{GMTVRB}} -K -O>> $outfile'.format(gmt_range, strsc), file=fout)
        print('if test \"$VHORIZONTAL\" -eq 1 ; then', file=fout)
        print('awk \'{{print $2,$3,$4,$5,$6,$7,$8,$1}}\' {:} | gmt psvelo -R -J -Se${{VSC}}/0.95/0 -W.3p,100 -A10p+e -V${{GMTVRB}} -Ggreen -O -K -L >> $outfile'.format(input_file), file=fout)
        print('fi', file=fout)
        print('if test \"$LABELS\" -eq 1 ; then', file=fout)
        print('awk \'{{print $2,$3,9,0,1,"RB",$1}}\' {:} | gmt pstext -Jm -R -Dj0.2c/0.2c -O -K -V${{GMTVRB}} >> $outfile'.format(input_file), file=fout)
        print('fi', file=fout)
        print('echo "9999 9999" | gmt psxy -J -R -O >> $outfile', file=fout)
    with open('.strain.info.dat', 'w') as fout:
        for idx, stn in enumerate(tensor_lst):
            dct = stn.info()
            clat, clon = utm2ell(stn.value_of('x'), stn.value_of('y'), utm_zone)
            # code lat lon Kmax sKmax Kmin sKmin Az sAz E sE gtot sgtot
            sigma = 1e-3
            print('{:} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}'.format(idx, degrees(clat), degrees(clon), dct['k_max'], sigma, dct['k_min'], sigma, dct['az'], sigma, dct['e'], sigma, dct['strain'], sigma), file=fout)
    with open('.station.info.dat', 'w') as fout:
        for idx, sta in enumerate(sta_lst):
            print('{:} {:} {:}'.format(sta.name, degrees(sta.lat), degrees(sta.lon)), file=fout)
    return

def plot_map(sta_list, stensor_list):
    lat0    = degrees(sum([ x.lat for x in sta_list ])/len(sta_list))
    lon0    = degrees(sum([ x.lon for x in sta_list ])/len(sta_list))
    lons    = [ degrees(x.lon) for x in sta_list ]
    lats    = [ degrees(x.lat) for x in sta_list ]
    lon_off = (max(lons)-min(lons))/10
    lat_off = (max(lats)-min(lats))/10
    my_map = Basemap(projection='merc', lat_0 = lat0, lon_0 = lon0, resolution = 'c', llcrnrlon=min(lons)-lon_off, llcrnrlat=min(lats)-lat_off, urcrnrlon=max(lons)+lon_off, urcrnrlat=max(lats)+lat_off)
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color = 'coral')
    my_map.drawmapboundary()
    my_map.drawmeridians(numpy.arange(floor(min(lons)), ceil(max(lons)), 2), labels=[True,False,False,True])
    my_map.drawparallels(numpy.arange(floor(min(lats)), ceil(max(lats)), 2), labels=[False,True,True,False], fontsize=10)

    for sta in sta_list:
        x, y = my_map(degrees(sta.lon), degrees(sta.lat))
        my_map.plot(x, y, 'bo', markersize=10)
        plt.text(x, y, sta.name)

    for tnr in stensor_list:
        x, y = my_map(degrees(tnr.lon), degrees(tnr.lat))
        my_map.plot(x, y, 'r+', markersize=8)

    print('[DEBUG] Area is {}/{}/{}/{}'.format(min(lons), max(lons), min(lats), max(lats)))
    plt.show()
    return

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description='Estimate Strain Tensor(s) from GNSS derived velocities.',
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
    required=True,
    help='The input file. This must be an ascii file containing the columns: \'station-name longtitude latitude Ve Vn SigmaVe SigmaVn Sne time-span\'. Longtitude and latitude must be given in decimal degrees; velocities (in east and north components) in mm/yr. Columns should be seperated by whitespaces. Note that at his point the last two columns (aka Sne and time-span) are not used, so they could have random values.')

parser.add_argument('--x-grid-step',
    default=50000,
    metavar='X_GRID_STEP',
    dest='x_grid_step',
    type=float,
    required=False,
    help='The x-axis grid step size in meters. This option is only relevant if the program computes more than one strain tensors.')

parser.add_argument('--y-grid-step',
    default=50000,
    metavar='Y_GRID_STEP',
    dest='y_grid_step',
    type=float,
    required=False,
    help='The y-axis grid step size in meters. This option is only relevant if the program computes more than one strain tensors.')

parser.add_argument('-m', '--method',
    default='shen',
    metavar='METHOD',
    dest='method',
    choices=['shen', 'veis'],
    required=False,
    help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')

parser.add_argument('-r', '--region',
    metavar='REGION',
    dest='region',
    help='Specify a region; any station (in the input file) falling outside will be ommited. The region should be given as a rectangle, specifying min/max values in longtitude and latitude (using decimal degrees). E.g. \"[...] --region=21.0/23.5/36.0/38.5 [...]\"',
    required=False)

parser.add_argument('-b', '--barycenter',
    dest='one_tensor',
    action='store_true',
    help='Only estimate one strain tensor, at the region\'s barycentre.')

parser.add_argument('--max-beta-angle',
    default=180,
    metavar='MAX_BETA_ANGLE',
    dest='max_beta_angle',
    type=float,
    required=False,
    help='Only relevant for \'--mehod=shen\'. Before estimating a tensor, the angles between consecutive points are computed. If the max angle is larger than max_beta_angle (in degrees), then the point is ommited (aka no tensor is computed). This option is used to exclude points from the computation tha only have limited geometric coverage (e.g. the edges of the grid).')

parser.add_argument('-t', '--weighting-function',
    default='gaussian',
    metavar='WEIGHTING_FUNCTION',
    dest='ltype',
    choices=['gaussian', 'quadratic'],
    required=False,
    help='Only relevant for \'--mehod=shen\'. Choose between a \'gaussian\' or a \'quadratic\' spatial weighting function.')

parser.add_argument('--Wt',
    default=24,
    metavar='Wt',
    dest='Wt',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. Let W=Σ_i*G_i, the total reweighting coefficients of the data, and let Wt be the threshold of W. For a given Wt, the smoothing constant D is determined by Wd=Wt . It should be noted that W is a function of the interpolation coordinate, therefore for the same Wt assigned, D varies spatially based on the in situ data strength; that is, the denser the local data array is, the smaller is D, and vice versa.')

parser.add_argument('--dmin',
    default=1,
    metavar='D_MIN',
    dest='dmin',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the lower limit for searching for an optimal d-param value. Unit is km.')

parser.add_argument('--dmax',
    default=500,
    metavar='D_MAX',
    dest='dmax',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the upper limit for searching for an optimal d-param value. Unit is km.')

parser.add_argument('--dstep',
    default=2,
    metavar='D_STEP',
    dest='dstep',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the step size for searching for an optimal d-param value. Unit is km.')

parser.add_argument('--min-spatial-weight',
    default=1e-6,
    metavar='MIN_SPATIAL_WEIGHT',
    dest='min_l_weight',
    type=float,
    required=False,
    help='Only relevant for \'--mehod=shen\'. Any station with a spatial weight smaller than this limit, will be excluded from the strain estimation (for a given point).')

parser.add_argument('--d-param',
    default=None,
    metavar='D_PARAMETER',
    dest='d_coef',
    type=float,
    required=False,
    help='Only relevant for \'--mehod=shen\'. This is the \'D\' parameter for computing the spatial weights. If this option is used, then the parameters: dmin, dmax, dstep and Wt are not used.')

##  Parse command line arguments.
args = parser.parse_args()
##  Collect args into a dictionary
dargs = vars(args)

##  Parse stations from input file; at input, station coordinates are in decimal
##+ degrees and velocities are in mm/yr.
##+ After reading, station coordinates are in radians and velocities are in
##+ m/yr.
sta_list_ell = parse_ascii_input(args.gps_file)
print('[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell)))

##  If a region is passed in, then only keep the stations that fall within it.
##  The region coordinates (min/max pairs) should be given in decimal degrees.
if args.region:
    try:
        Napr = len(sta_list_ell)
        lonmin, lonmax, latmin, latmax = [ radians(float(i)) for i in args.region.split('/') ]
        sta_list_ell = cut_rectangle(lonmin, lonmax, latmin, latmax, sta_list_ell)
        Npst = len(sta_list_ell)
        print('[DEBUG] Station filtered to fit input region: {:7.3f}/{:7.3f}/{:7.3f}/{:7.3f}'.format(lonmin, lonmax, latmin, latmax))
        print('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))
    except:
        ## TODO we should exit with error here
        print('[ERROR] Failed to parse region argument \"{:}\"'.format(args.region))

##  Make a new station list (copy of the original one), where all coordinates
##+ are in UTM. All points should belong to the same ZONE.
##  Note that station ellipsoidal coordinates are in radians while the cartesian
##+ coordinates are in meters.
mean_lon = degrees(sum([ x.lon for x in sta_list_ell ]) / len(sta_list_ell))
utm_zone = floor(mean_lon/6)+31
utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
print('[DEBUG] Mean longtitude is {} deg.; using Zone = {} for UTM'.format(mean_lon, utm_zone))
sta_list_utm = deepcopy(sta_list_ell)
for idx, sta in enumerate(sta_list_utm):
    N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), utm_zone)
    sta_list_utm[idx].lon = E
    sta_list_utm[idx].lat = N
    assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
print('[DEBUG] Station list transformed to UTM.')

##  Compute only one Strain Tensor, at the region's barycenter; then exit.
if args.one_tensor:
    if args.method == 'shen':
        sstr = ShenStrain(0e0, 0e0, sta_list_utm, **dargs)
    else:
        sstr = VeisStrain(0e0, 0e0, sta_list_utm)
    sstr.set_to_barycenter()
    sstr.estimate()
    #gmt_script(sta_list_ell, [sstr], utm_zone, outfile='gmt_script', projscale=2000000, strsc=5, frame=2)
    sys.exit(0)

##  Construct the grid, based on station coordinates (Ref. UTM)
fout = open('strain_info.dat', 'w')
print('{:^9s} {:^9s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s}'.format('X (m)', 'Y (m)', 'Ux', 'Uy', 'tauy', 'omega', 'taux', 'tauxy', 'emax', 'emin', 'taumax', 'dexazim', 'dilat'), file=fout)
strain_list = []
if args.method == 'shen':
    grd = pystrain.grid.generate_grid(sta_list_utm, args.x_grid_step, args.y_grid_step)
    print('[DEBUG] Constructed the grid. Limits are:')
    print('\tEasting : from {} to {} with step {}'.format(grd.x_min, grd.x_max, grd.x_step))
    print('\tNorthing: from {} to {} with step {}'.format(grd.y_min, grd.y_max, grd.y_step))
    print('[DEBUG] Estimating strain tensor for each cell center')
    ##  Iterate through the grid (on each cell center)
    node_nr = 0
    for x, y in grd:
        clat, clon = utm2ell(x, y, utm_zone)
        print('[DEBUG] Grid point at {:7.4f}, {:7.4f} or {:}, {:}'.format(degrees(clon), degrees(clat), x, y))
        sstr = ShenStrain(x, y, sta_list_utm, **dargs)
        if degrees(max(sstr.beta_angles())) <= args.max_beta_angle:
            try:
                estim2 = sstr.estimate()
                node_nr += 1
                print('[DEBUG] Computed tensor at {:7.4f}, {:7.4f} for node {:3d}/{:3d}'.format(degrees(clon), degrees(clat), node_nr, grd.xpts*grd.ypts))
                sstr.print_details(fout)
                strain_list.append(sstr)
            except RuntimeError:
                print('[DEBUG] Too few observations to estimate strain at {:7.4f}, {:7.4f}'.format(degrees(clon), degrees(clat)))
        else:
            print('[DEBUG] Skipping computation at {:7.4f},{:7.4f} because of limited coverage (max_beta= {:6.2f}deg.)'.format(degrees(clon), degrees(clat), degrees(max(sstr.beta_angles()))))
else:
    points = numpy.array([ [sta.lon, sta.lat] for sta in sta_list_utm ])
    tri = Delaunay(points)
    for idx, trng in enumerate(tri.simplices):
        cx = (sta_list_utm[trng[0]].lon + sta_list_utm[trng[1]].lon + sta_list_utm[trng[2]].lon)/3e0
        cy = (sta_list_utm[trng[0]].lat + sta_list_utm[trng[1]].lat + sta_list_utm[trng[2]].lat)/3e0
        sstr = VeisStrain(cx, cy, [sta_list_utm[trng[0]], sta_list_utm[trng[1]], sta_list_utm[trng[2]]])
        estim2 = sstr.estimate()
        clat, clon = utm2ell(cx, cy, utm_zone)
        strain_list.append(sstr)

fout.close()
gmt_script(args.gps_file, sta_list_ell, strain_list, utm_zone, outfile='gmt_script', projscale=2000000, strsc=5, frame=2)
#plot_map(sta_list_ell, strain_list)
