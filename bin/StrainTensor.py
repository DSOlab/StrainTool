#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
############################################## standard libs
import sys
import os
import time
from datetime import datetime
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
from scipy.spatial import Delaunay

Version = 'StrainTensor.py Version: 1.0-beta (pre-release)'

def cut_rectangle(xmin, xmax, ymin, ymax, sta_lst, sta_list_to_degrees=False):
    new_sta_lst = []
    for sta in sta_lst:
        if sta_list_to_degrees:
            slon = degrees(sta.lon)
            slat = degrees(sta.lat)
        else:
            slon = sta.lon
            slat = sta.lat
        if slon >= xmin and slon <= xmax and slat >= ymin and slat <= ymax:
            new_sta_lst.append(sta)
    return new_sta_lst

def write_station_info(sta_lst, filename='station_info.dat'):
    with open(filename, 'w') as fout:
        print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format('Station', 'Longtitude', 'Latitude', 'Ve', 'Vn', 'sVe', 'sVn'), file=fout)
        print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format('', 'deg.', 'deg', 'mm/yr', 'mm/yr', 'mm/yr', 'mm/yr'), file=fout)
        for idx, sta in enumerate(sta_lst):
            print('{:10s} {:+10.5f} {:10.5f} {:+7.2f} {:+7.2f} {:+7.3f} {:+7.3f}'.format(sta.name, degrees(sta.lon), degrees(sta.lat), sta.ve*1e03, sta.vn*1e03, sta.se*1e03, sta.sn*1e03), file=fout)
    return

def print_model_info(fout, cmd, clargs):
    print('{:}'.format(Version), file=fout)
    print('Command used:\n\t{:}'.format(' '.join(cmd)), file=fout)
    print('Run at: {:}'.format(datetime.now().strftime('%c')), file=fout)
    print('Command line switches/options parsed:', file=fout)
    for key in clargs:
        print('\t{:20s} -> {:}'.format(key, clargs[key]), file=fout)
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
    default=0.5,
    metavar='X_GRID_STEP',
    dest='x_grid_step',
    type=float,
    required=False,
    help='The x-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors. Default is 0.5(deg).')

parser.add_argument('--y-grid-step',
    default=0.5,
    metavar='Y_GRID_STEP',
    dest='y_grid_step',
    type=float,
    required=False,
    help='The y-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors. Default is 0.5(deg)')

parser.add_argument('-m', '--method',
    default='shen',
    metavar='METHOD',
    dest='method',
    choices=['shen', 'veis'],
    required=False,
    help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter. Default is \'shen\'.')

parser.add_argument('-r', '--region',
    metavar='REGION',
    dest='region',
    help='Specify a region; any station (in the input file) falling outside will be ommited. The region should be given as a rectangle, specifying min/max values in longtitude and latitude (using decimal degrees). E.g. \"[...] --region=21.0/23.5/36.0/38.5 [...]\"',
    required=False)

parser.add_argument('-c', '--cut-excess-stations',
    dest='cut_outoflim_sta',
    help='This option is only considered if the \'-r\' option is set. If this this option is enabled, then any station (from the input file) outside the region limit (passed in via the \'-r\' option) is not considered in the strain estimation.',
    action='store_true')

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
    help='Only relevant for \'--mehod=shen\'. Before estimating a tensor, the angles between consecutive points are computed. If the max angle is larger than max_beta_angle (in degrees), then the point is ommited (aka no tensor is computed). This option is used to exclude points from the computation tha only have limited geometric coverage (e.g. the edges of the grid). Default is 180 deg.')

parser.add_argument('-t', '--weighting-function',
    default='gaussian',
    metavar='WEIGHTING_FUNCTION',
    dest='ltype',
    choices=['gaussian', 'quadratic'],
    required=False,
    help='Only relevant for \'--mehod=shen\'. Choose between a \'gaussian\' or a \'quadratic\' spatial weighting function. Default is \'gaussian\'.')

parser.add_argument('--Wt',
    default=24,
    metavar='Wt',
    dest='Wt',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. Let W=Σ_i*G_i, the total reweighting coefficients of the data, and let Wt be the threshold of W. For a given Wt, the smoothing constant D is determined by Wd=Wt . It should be noted that W is a function of the interpolation coordinate, therefore for the same Wt assigned, D varies spatially based on the in situ data strength; that is, the denser the local data array is, the smaller is D, and vice versa. Default is Wt=24.')

parser.add_argument('--dmin',
    default=1,
    metavar='D_MIN',
    dest='dmin',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the lower limit for searching for an optimal d-param value. Unit is km. Default is dmin=1km.')

parser.add_argument('--dmax',
    default=500,
    metavar='D_MAX',
    dest='dmax',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the upper limit for searching for an optimal d-param value. Unit is km. Default is dmax=500km.')

parser.add_argument('--dstep',
    default=2,
    metavar='D_STEP',
    dest='dstep',
    type=int,
    required=False,
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the step size for searching for an optimal d-param value. Unit is km. Default is dstep=2km.')

parser.add_argument('--d-param',
    default=None,
    metavar='D_PARAMETER',
    dest='d_coef',
    type=float,
    required=False,
    help='Only relevant for \'--mehod=shen\'. This is the \'D\' parameter for computing the spatial weights. If this option is used, then the parameters: dmin, dmax, dstep and Wt are not used.')

parser.add_argument('-g', '--generate-statistics',
    dest='generate_stats',
    help='Only relevant when \'--mehod=shen\' and \'--barycenter\' is not set. This option will create an output file, named \'strain_stats.dat\', where estimation info and statistics will be written.',
    action='store_true')

parser.add_argument('--verbose',
    dest='verbose_mode',
    help='Run in verbose mode (show debugging messages)',
    action='store_true')

parser.add_argument('-v',
    dest='version',
    help='Display version and exit.',
    action='store_true')

##  Time the program (for opt/ing purpose only)
start_time = time.time()

##  Parse command line arguments and stack them in a dictionary
args  = parser.parse_args()
dargs = vars(args)

##  Wait!! maybe we only want the version
if args.version:
    print('{}'.format(Version))
    sys.exit(0)

## Verbose print (function only exists in verbose mode)
vprint = print if args.verbose_mode else lambda *a, **k: None

## If needed, open a file to write model info and statistics
fstats = open('strain_stats.dat', 'w') if args.generate_stats else None
if fstats: print_model_info(fstats, sys.argv, dargs)

##  Parse stations from input file; at input, station coordinates are in decimal
##+ degrees and velocities are in mm/yr.
##  After reading, station coordinates are in radians and velocities are in
##+ m/yr.
if not os.path.isfile(args.gps_file):
    print('[ERROR] Cannot find input file \'{}\'.'.format(args.gps_file), file=sys.stderr)
    sys.exit(1)
sta_list_ell = parse_ascii_input(args.gps_file)
print('[DEBUG] Reading station coordinates and velocities from {}'.format(args.gps_file))
print('[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell)))

##  If a region is passed in, resolve it.
##+ If cutting out-of-limits stations option is set, or method is veis, then 
##+ only keep the stations that fall within it.
##  The region coordinates (min/max pairs) should be given in decimal degrees.
if args.region:
    try:
        lonmin, lonmax, latmin, latmax = [ float(i) for i in args.region.split('/') ]
        if args.cut_outoflim_sta or args.method == 'veis':
            Napr = len(sta_list_ell)
            #  Note that we have to convert radians to degrees for station 
            #+ coordinates, hence 'sta_list_to_degrees=True'
            sta_list_ell = cut_rectangle(lonmin, lonmax, latmin, latmax, sta_list_ell, True)
            Npst = len(sta_list_ell)
            vprint('[DEBUG] Stations filtered to fit input region: {:7.3f}/{:7.3f}/{:7.3f}/{:7.3f}'.format(lonmin, lonmax, latmin, latmax))
            vprint('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))
    except:
        ## TODO we should exit with error here
        print('[ERROR] Failed to parse region argument \"{:}\"'.format(args.region), file=sys.stderr)

##  Filter out stations that are never going to be used. This is an opt!
if args.region and not args.method == 'veis' and not args.cut_outoflim_sta:
    vprint('[DEBUG] Filtering stations based on their distance from region barycentre.')
    Napr = len(sta_list_ell)
    mean_lon, mean_lat = radians(lonmin+(lonmax-lonmin)/2e0), radians(latmin+(latmax-latmin)/2e0)
    # print('[DEBUG] Barycentre set to {:+10.5f}, {:10.5f}'.format(lonmin+(lonmax-lonmin)/2e0, latmin+(latmax-latmin)/2e0))
    bc =  Station(lon=mean_lon, lat=mean_lat)
    endpt = Station(lon=radians(lonmax), lat=radians(latmax))
    cutoffdis = abs(endpt.haversine_distance(bc)/1e3) # barycentre to endpoint (km)
    # print('[DEBUG] Distance from barycentre to endpoint is {:10.3f}km'.format(cutoffdis))
    d = 2e0*(args.d_coef if args.d_coef is not None else args.dmax)
    cutoffdis += d * (2.15e0 if args.ltype == 'gaussian' else 10e0) # in km
    vprint('[DEBUG] Using cut-off distance {:10.3f}km'.format(cutoffdis))
    #for s in sta_list_ell:
    #    d = s.haversine_distance(bc)
    #    print('station {:} is {:7.1f}km away {:}'.format(s.name, d/1e3, 'Acc' if s.haversine_distance(bc)/1e3 <= cutoffdis else 'Rej'))
    sta_list_ell = [ s for s in sta_list_ell if s.haversine_distance(bc)/1e3 <= cutoffdis ]
    Npst = len(sta_list_ell)
    print('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))
    Npst = len(sta_list_ell)

##  Make a new station list (copy of the original one), where all coordinates
##+ are in UTM. All points should belong to the same ZONE.
##  Note that station ellipsoidal coordinates are in radians while the cartesian
##+ coordinates are in meters.
##
##  TODO is this mean_lon the optimal?? or should it be the region's mean longtitude
##
mean_lon = degrees(sum([ x.lon for x in sta_list_ell ]) / len(sta_list_ell))
utm_zone = floor(mean_lon/6)+31
utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
vprint('[DEBUG] Mean longtitude is {} deg.; using Zone = {} for UTM'.format(mean_lon, utm_zone))
sta_list_utm = deepcopy(sta_list_ell)
for idx, sta in enumerate(sta_list_utm):
    N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), utm_zone)
    sta_list_utm[idx].lon = E
    sta_list_utm[idx].lat = N
    assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
vprint('[DEBUG] Station list transformed to UTM.')

##  Open file to write Strain Tensor estimates; write the header
fout = open('strain_info.dat', 'w')
vprint('[DEBUG] Strain info written in file: {}'.format('strain_info.dat'))
print('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}'.format('Latitude', 'Longtitude', 'vx+dvx', 'vy+dvy', 'w+dw', 'exx+dexx', 'exy+dexy', 'eyy+deyy', 'emax+demax', 'emin+demin', 'shr+dshr', 'azi+dazi', 'dilat+ddilat', 'sec. invariant'), file=fout)
print('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}'.format('deg', 'deg', 'mm/yr', 'mm/yr', 'nrad/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'deg.', 'nstrain/yr', 'nstrain/yr'), file=fout)

##  Compute only one Strain Tensor, at the region's barycenter; then exit.
if args.one_tensor:
    print('[DEBUG] Estimating Strain Tensor at region\'s barycentre.')
    if args.method == 'shen':
        sstr = ShenStrain(0e0, 0e0, sta_list_utm, **dargs)
    sstr.set_to_barycenter()
    sstr.estimate()
    sys.exit(0)

# strain_list = [] Probably we do not need to keep the tensors ...
if args.method == 'shen':  ## Going for Shen algorithm ...
    ##  Construct the grid, in ellipsoidal coordinates --degrees--. If a region
    ##+ is not passed in, the grid.generate_grid will transform lon/lat pairs 
    ##+ to degrees and produce a grid from extracting min/max crds from the
    ##+ station list.
    if args.region:
        grd = pystrain.grid.Grid(lonmin, lonmax, args.x_grid_step, latmin, latmax, args.y_grid_step)
    else:
        grd = pystrain.grid.generate_grid(sta_list_ell, args.x_grid_step, args.y_grid_step, True)
    print('[DEBUG] Grid Information:')
    print('[DEBUG]\tLongtitude : from {} to {} with step {} (deg)'.format(grd.x_min, grd.x_max, grd.x_step))
    print('[DEBUG]\tLatitude   : from {} to {} with step {} (deg)'.format(grd.y_min, grd.y_max, grd.y_step))
    print('[DEBUG] Number of Strain Tensors to be estimated: {}'.format(grd.xpts*grd.ypts))
    if fstats: print('{:8s} {:8s} {:12s} {:12s} {:12s} {:7s}'.format('Longtitude','Latitude','# stations', 'D (optimal)','CutOff dis.', 'Sigma'), file=fstats)
    vprint('[DEBUG] Estimating strain tensor for each cell center:')
    ##  Iterate through the grid (on each cell center). Grid returns cell-centre
    ##+ coordinates in lon/lat pairs, in degrees!
    node_nr, nodes_estim = 0, 0
    for x, y in grd:
        clat, clon =  radians(y), radians(x)
        N, E, ZN, _ = ell2utm(clat, clon, Ellipsoid("wgs84"), utm_zone)
        assert ZN == utm_zone
        vprint('[DEBUG] Grid point at {:+8.4f}, {:8.4f} or E={:}, N={:}'.format(x, y, E, N))
        print('[DEBUG] {:5d}/{:7d}'.format(node_nr+1, grd.xpts*grd.ypts), end="\r")
        ## Construct the Strain instance, with all args (from input)
        sstr = ShenStrain(E, N, sta_list_utm, **dargs)
        ## check azimouth coverage (aka max β angle)
        if degrees(max(sstr.beta_angles())) <= args.max_beta_angle:
            try:
                sstr.estimate()
                vprint('[DEBUG] Computed tensor at {:+8.4f}, {:8.4f} for node {:3d}/{:3d}'.format(x, y, node_nr+1, grd.xpts*grd.ypts))
                sstr.print_details(fout, utm_zone)
                if fstats: print('{:+8.4f}, {:8.4f} {:12d} {:12.2f} {:12.2f} {:7.3f}'.format(x,y,len(sstr.__stalst__), sstr.__options__['d_coef'],sstr.__options__['cutoff_dis'], sstr.__sigma0__), file=fstats)
                # strain_list.append(sstr)
                nodes_estim += 1
            except RuntimeError:
                vprint('[DEBUG] Too few observations to estimate strain at {:+8.4f}, {:8.4f}. Point skipped.'.format(x,y))
            except ArithmeticError:
                vprint('[DEBUG] Failed to compute parameter VcV matrix for strain at {:+8.4f}, {:8.4f}. Point skipped'.format(x,y))
        else:
            vprint('[DEBUG] Skipping computation at {:+8.4f},{:8.4f} because of limited coverage (max_beta= {:6.2f}deg.)'.format(x, y, degrees(max(sstr.beta_angles()))))
        node_nr += 1
    print('[DEBUG] Estimated Strain Tensors for {} out of {} nodes'.format(nodes_estim, node_nr))
else:
    ## Open file to write delaunay triangles.
    print('[DEBUG] Estimating Strain Tensors at the barycentre of Delaunay triangles')
    dlnout = open('delaunay_info.dat', 'w')
    points = numpy.array([ [sta.lon, sta.lat] for sta in sta_list_utm ])
    tri = Delaunay(points)
    print('[DEBUG] Number of Delaunay triangles: {}'.format(len(tri.simplices)))
    for idx, trng in enumerate(tri.simplices):
        print('[DEBUG] {:5d}/{:7d}'.format(idx+1, len(tri.simplices)), end="\r")
        ## triangle barycentre
        cx = (sta_list_utm[trng[0]].lon + sta_list_utm[trng[1]].lon + sta_list_utm[trng[2]].lon)/3e0
        cy = (sta_list_utm[trng[0]].lat + sta_list_utm[trng[1]].lat + sta_list_utm[trng[2]].lat)/3e0
        ##  Construct a strain instance, at the triangle's barycentre, with only
        ##+ 3 points (in UTM) and equal_weights weighting scheme.
        sstr = ShenStrain(cx, cy, [sta_list_utm[trng[0]], sta_list_utm[trng[1]], sta_list_utm[trng[2]]], weighting_function='equal_weights')
        sstr.estimate()
        sstr.print_details(fout, utm_zone)
        ## Print the triangle in the corresponding file (ellipsoidal crd, degrees)
        print('> {:}, {:}, {:}'.format(sta_list_utm[trng[0]].name, sta_list_utm[trng[1]].name, sta_list_utm[trng[2]].name), file=dlnout)
        print('{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}'.format(*[ degrees(x) for x in [sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat, sta_list_ell[trng[1]].lon, sta_list_ell[trng[1]].lat, sta_list_ell[trng[2]].lon, sta_list_ell[trng[2]].lat, sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat]]), file=dlnout)
        # strain_list.append(sstr)
    dlnout.close()

fout.close()
write_station_info(sta_list_ell)
print('[DEBUG] Total running time: {:10.2f} sec.'.format((time.time() - start_time)))
