#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import time
from datetime import datetime
from copy import deepcopy
from math import degrees, radians, floor, ceil
import numpy
from scipy.spatial import Delaunay
import argparse
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *
import pystrain.grid

Version = 'StrainTensor.py Version: 1.0-r1'
STRAIN_OUT_FILE = 'strain_info.dat'
STATISTICS_FILE = 'strain_stats.dat'

def cut_rectangle(xmin, xmax, ymin, ymax, sta_lst, sta_list_to_degrees=False):
    """ Filter stations that are located within a rectange. The rectangle is
        +-------+--ymax
        |       |
        +-------+--ymin
        |       |
        xmin    xmax 
        The function will return a new list, where for each of the stations,
        the following is true:
            * xmin <= station.lon <= xmax and
            * ymin <= station.lat <= ymax
        If the argument 'sta_list_to_degrees' is set to True, then before
        comaring, each of the station's lon and lat are transformed to degrees
        (they are supposed to be in radians).
    """
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
    """ Write station information to an output file. sta_list if a list of
        Stations.
        Station information are written as:
        Station Longtitude Latitude Ve Vn sVe sVn
                deg.       deg      mm/yr
        The file to be written is named as $filename
    """
    with open(filename, 'w') as fout:
        print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format(
            'Station', 'Longtitude', 'Latitude', 'Ve', 'Vn', 'sVe', 'sVn'),
            file=fout)
        print('{:^10s} {:^10s} {:^10s} {:7s} {:7s} {:7s} {:7s}'.format(
            '', 'deg.', 'deg', 'mm/yr', 'mm/yr', 'mm/yr', 'mm/yr'),
            file=fout)
        for idx, sta in enumerate(sta_lst):
            print('{:10s} {:+10.5f} {:10.5f} {:+7.2f} {:+7.2f} {:+7.3f} {:+7.3f}'.format(
                sta.name, degrees(sta.lon), degrees(sta.lat), sta.ve*1e03, 
                sta.vn*1e03, sta.se*1e03, sta.sn*1e03), file=fout)
    return

def print_model_info(fout, cmd, clargs):
    """ Write basic information to an open output stream (e.g. a file).
    """
    print('{:}'.format(Version), file=fout)
    print('Command used:\n\t{:}'.format(' '.join(cmd)), file=fout)
    print('Run at: {:}'.format(datetime.now().strftime('%c')), file=fout)
    print('Command line switches/options parsed:', file=fout)
    for key in clargs:
        print('\t{:20s} -> {:}'.format(key, clargs[key]), file=fout)
    return

def compute__(igrd, sta_list_utm, utm_lcm, fout, fstats, vprint_fun, **dargs):
    """ Function to perform the bulk of a Strain Tensor estimation.
        For each of the grid cells, a ShenStrain object will be created, using
        the list of stations and the **dargs options.

        Args:
            grd (pystrain::Grid): The grid; one straintensor per cell is
                                  estimated (at the centre of the grid)
            sta_list_utm (list of Station): The list of stations to be used for
                                  strain tensor estimation
            utmzone (float):      The UTM zone used to convert ellipsoidal to
                                  UTM coordinates.
            fout (output stream): An (open) output stream where estimation results
                                  (aka strain information) are to be written
            fstats (output stream): An (open) output stream where estimation
                                  statistics are written
            vprint_fun (function) : A function that handles printing. Based on
                                  user options we may want or not to print
                                  verbose information. This function does exactly
                                  that. Normally, this function is just the
                                  normal print function or a no-opt, see vprint(...)
                                  defined in __main__
            **dargs (dictionary)  : A list of parameters to use when constructing
                                  the individual Strain Tensors

        Warning:
            The output streams are passed in open but are closed by the function!
            Leaving the streams open, may cause not proper reporting of results
            in Python v2.x and in multithreading mode (probably the streams are 
            not flushed before returning or something). Anyway, always close the 
            streams before exiting.
    """
    #print('--> Thread given grid : X:{:}/{:}/{:} Y:{:}/{:}/{:}'.format(igrd.x_min, igrd.x_max, igrd.x_step, igrd.y_min, igrd.y_max, igrd.y_step))
    node_nr, nodes_estim = 0, 0
    for x, y in igrd:
        clat, clon =  radians(y), radians(x)
        #print('--> computing tensor at lon {:}, lat {:}'.format(x, y))
        N, E, ZN, lcm = ell2utm(clat, clon, Ellipsoid("wgs84"), utm_lcm)
        #assert ZN == utmzone
        assert utm_lcm == lcm
        vprint_fun('[DEBUG] Grid point at {:+8.4f}, {:8.4f} or E={:}, N={:}'.format(
            x, y, E, N))
        if not dargs['multiproc_mode']:
            print('[DEBUG] {:5d}/{:7d}'.format(node_nr+1, grd.xpts*grd.ypts), end="\r")
        ## Construct the Strain instance, with all args (from input)
        # sstr = ShenStrain(E, N, sta_list_utm, **dargs)
        sstr = ShenStrain(E, N, clat<0e0, sta_list_utm, **dargs)
        ## check azimouth coverage (aka max β angle)
        if degrees(max(sstr.beta_angles())) <= dargs['max_beta_angle']:
            try:
                sstr.estimate()
                vprint_fun('[DEBUG] Computed tensor at {:+8.4f} {:+8.4f} for node {:3d}/{:3d}'.format(x, y, node_nr+1, grd.xpts*grd.ypts))
                sstr.print_details_v2(fout, utm_lcm)
                if fstats: print('{:+9.4f} {:+10.4f} {:6d} {:14.2f} {:10.2f} {:12.3f}'.format(x,y,len(sstr.__stalst__), sstr.__options__['d_coef'],sstr.__options__['cutoff_dis'], sstr.__sigma0__), file=fstats)
                nodes_estim += 1
            except RuntimeError:
                vprint_fun('[DEBUG] Too few observations to estimate strain at {:+8.4f}, {:8.4f}. Point skipped.'.format(x,y))
            except ArithmeticError:
                vprint_fun('[DEBUG] Failed to compute parameter VcV matrix for strain at {:+8.4f}, {:8.4f}. Point skipped'.format(x,y))
        else:
            vprint_fun('[DEBUG] Skipping computation at {:+8.4f},{:8.4f} because of limited coverage (max_beta= {:6.2f}deg.)'.format(x, y, degrees(max(sstr.beta_angles()))))
        node_nr += 1
    print('[DEBUG] Estimated Strain Tensors for {} out of {} nodes'.format(nodes_estim, node_nr))
    fout.close()
    if fstats: fstats.close()

##  If only the formatter_class could be:
##+ argparse.RawTextHelpFormatter|ArgumentDefaultsHelpFormatter ....
##  Seems to work with multiple inheritance!
class myFormatter(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description='Estimate Strain Tensor(s) from GNSS derived velocities.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
Send bug reports to:
  Xanthos Papanikolaou, xanthos@mail.ntua.gr
  Dimitris Anastasiou,dganastasiou@gmail.com
September, 2021'''))

parser.add_argument('-i', '--input-file',
    default=argparse.SUPPRESS,
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
    help='The x-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')

parser.add_argument('--y-grid-step',
    default=0.5,
    metavar='Y_GRID_STEP',
    dest='y_grid_step',
    type=float,
    required=False,
    help='The y-axis grid step size in degrees. This option is only relevant if the program computes more than one strain tensors.')

parser.add_argument('-m', '--method',
    default='shen',
    metavar='METHOD',
    dest='method',
    choices=['shen', 'veis'],
    required=False,
    help='Choose a method for strain estimation. If \'shen\' is passed in, the estimation will follow the algorithm described in Shen et al, 2015, using a weighted least squares approach. If \'veis\' is passed in, then the region is going to be split into delaneuy triangles and a strain estimated in each barycenter.')

parser.add_argument('-r', '--region',
    default=argparse.SUPPRESS,
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
    help='Only relevant for \'--mehod=shen\' and if \'d-param\' is not passed in. This is the lower limit for searching for an optimal D-parameter value. Unit is km.')

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

parser.add_argument('--multicore',
    dest='multiproc_mode',
    help='Run in multithreading mode',
    action='store_true')

parser.add_argument('-v',
    dest='version',
    help='Display version and exit.',
    action='store_true')

if __name__ == '__main__':
    ##  Wait!! maybe the user just paseed in "-v" without an input file. Do not
    ##+ resolve the parser yet (it ll cause an error)
    if len(sys.argv[1:]) == 1 and sys.argv[1] == "-v":
        print('{}'.format(Version))
        sys.exit(0)

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

    ## if in mutlithreading mode, load the module
    if args.multiproc_mode:
        if args.method == 'shen':
            import multiprocessing
            cpu_count = multiprocessing.cpu_count()
            print("[DEBUG] Using multithreaded version; available CPU's: {:02d}".format(
                cpu_count))
        else:
            print("[DEBUG] Multithreading is only available when using shen method; ignoring the \"--multicore\" switch!")
    
	## import dill module for windows multithreading processing
    if args.multiproc_mode and os.name == 'nt':
        print("[DEBUG] Import dill module for windows multithreading processing")
        import dill
        

    ## If needed, open a file to write model info and statistics
    fstats = open(STATISTICS_FILE, 'w') if args.generate_stats else None
    if fstats: print_model_info(fstats, sys.argv, dargs)

    ##  Parse stations from input file; at input, station coordinates are in decimal
    ##+ degrees and velocities are in mm/yr.
    ##  After reading, station coordinates are in radians and velocities are in
    ##+ m/yr.
    if not os.path.isfile(args.gps_file):
        print('[ERROR] Cannot find input file \'{}\'.'.format(
            args.gps_file), file=sys.stderr)
        sys.exit(1)
    try:
        sta_list_ell = parse_ascii_input(args.gps_file, args.method=='shen')
    except ValueError as err:
        print(err)
        print('[ERROR] Failed to parse input file: \"{:}\"'.format(args.gps_file))
        sys.exit(1)
    print('[DEBUG] Reading station coordinates and velocities from {}'.format(
        args.gps_file))
    print('[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell)))

    ##  If a region is passed in, resolve it (from something like 
    ##+ '21.0/23.5/36.0/38.5'). Note that limits are in dec. degrees.
    ##+ If cutting out-of-limits stations option is set, or method is veis, then 
    ##+ only keep the stations that fall within it.
    ##  The region coordinates (min/max pairs) should be given in decimal degrees.
    if 'region' in args:
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
                if Npst < 3:
                    print('[DEBUG] Left with only {:d} stations! Cannot do anything'.format(Npst))
                    sys.exit(0)
        except:
            ## TODO we should exit with error here
            print('[ERROR] Failed to parse region argument \"{:}\"'.format(
                args.region), file=sys.stderr)

    ##  Filter out stations that are never going to be used. This is an opt!
    ##  This is only needed when the used has specified:
    ##+ '[...] --region=a/b/c/d --method='shen' [...]' and NOT --cut-excess-station
    ##+ because:
    ##+ * If there is no region, we suppose that we want all the region covered
    ##+   by the stations
    ##+ * If method='veis' we are using Delaneuey triangles anyway
    ##+ * If '--cut-excess-station' is set, we have already cut-off any stations
    ##+   outside the wanted region
    ##  This is performed as follows:
    ##+ 1. Compute distance from centre of region to point (lonmax, latmax), aka R
    ##+ 2. Compute D: User has specified 'D_PARAMETER'? D=2*D_PARAMETER else D=2*D_MAX
    ##+ 3. Compute C: WEIGHTING_FUNCTION='gaussian'? C=R+D*2.15 else C=R+D*10
    ##+ 4. Filter out any station that has distance from the centre > C
    ##  Note that all distances are computed via the Haversine formula and all units
    ##+ are Km
    if 'region' in args and not args.method == 'veis' and not args.cut_outoflim_sta:
        vprint('[DEBUG] Filtering stations based on their distance from region barycentre.')
        Napr = len(sta_list_ell)
        mean_lon, mean_lat = radians(lonmin+(lonmax-lonmin)/2e0), radians(latmin+(latmax-latmin)/2e0)
        bc =  Station(lon=mean_lon, lat=mean_lat)
        endpt = Station(lon=radians(lonmax), lat=radians(latmax))
        cutoffdis = abs(endpt.haversine_distance(bc)/1e3) # barycentre to endpoint (km)
        d = 2e0*(args.d_coef if args.d_coef is not None else args.dmax)
        cutoffdis += d * (2.15e0 if args.ltype == 'gaussian' else 10e0) # in km
        vprint('[DEBUG] Using cut-off distance {:10.3f}km'.format(cutoffdis))
        sta_list_ell = [ s for s in sta_list_ell if s.haversine_distance(bc)/1e3 <= cutoffdis ]
        Npst = len(sta_list_ell)
        print('[DEBUG] {:4d} out of original {:4d} stations remain to be processed.'.format(Npst, Napr))

    ##  Make a new station list (copy of the original one), where all coordinates
    ##+ are in UTM. All points should belong to the same ZONE.
    ##  Note that station ellipsoidal coordinates are in radians while the 
    ##+ cartesian (projection) coordinates are in meters.
    ##
    ##  TODO is this mean_lon the optimal?? or should it be the region's mean longtitude
    ##
    mean_lon = degrees(sum([ x.lon for x in sta_list_ell ]) / len(sta_list_ell))
    #utm_zone = floor(mean_lon/6)+31
    #utm_zone = utm_zone + int(utm_zone<=0)*60 - int(utm_zone>60)*60
    lcm = radians(floor(mean_lon))
    #print('[DEBUG] Mean longtitude is {} deg.; using Zone = {} for UTM'.format(mean_lon, utm_zone))
    sta_list_utm = deepcopy(sta_list_ell)
    for idx, sta in enumerate(sta_list_utm):
        N, E, Zone, lcm = ell2utm(sta.lat, sta.lon, Ellipsoid("wgs84"), lcm)
        sta_list_utm[idx].lon = E
        sta_list_utm[idx].lat = N
        # assert Zone == utm_zone, "[ERROR] Invalid UTM Zone."
    vprint('[DEBUG] Station list transformed to UTM.')

    ##  Open file to write Strain Tensor estimates; write the header
    fout = open(STRAIN_OUT_FILE, 'w')
    vprint('[DEBUG] Strain info written in file: {}'.format(STRAIN_OUT_FILE))
    print('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}'.format('Latitude', 'Longtitude', 'vx+dvx', 'vy+dvy', 'w+dw', 'exx+dexx', 'exy+dexy', 'eyy+deyy', 'emax+demax', 'emin+demin', 'shr+dshr', 'azi+dazi', 'dilat+ddilat', 'sec. invariant+dsec inv.'), file=fout)
    print('{:^9s} {:^9s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s} {:^15s}'.format('deg', 'deg', 'mm/yr', 'mm/yr', 'deg/Myr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'deg.', 'nstrain/yr', 'nstrain/yr'), file=fout)

    ##  Compute only one Strain Tensor, at the region's barycenter; then exit.
    if args.one_tensor:
        print('[DEBUG] Estimating Strain Tensor at region\'s barycentre.')
        if args.method == 'shen':
            sstr = ShenStrain(0e0, 0e0, False, sta_list_utm, **dargs)
        else:
            sstr = ShenStrain(0e0, 0e0, False, sta_list_utm, weighting_function='equal_weights')
        sstr.set_to_barycenter()
        sstr.estimate()
        sstr.print_details(fout, utm_lcm)
        fout.close()
        write_station_info(sta_list_ell)
        print('[DEBUG] Total running time: {:10.2f} sec.'.format((time.time() - start_time)))      
        sys.exit(0)

    if args.method == 'shen':  ## Going for Shen algorithm ...
        ##  Construct the grid, in ellipsoidal coordinates --degrees--. If a region
        ##+ is not passed in, the grid.generate_grid will transform lon/lat pairs 
        ##+ to degrees and produce a grid from extracting min/max crds from the
        ##+ station list.
        if 'region' in args:
            grd = pystrain.grid.Grid(lonmin, lonmax, args.x_grid_step, latmin, latmax, args.y_grid_step)
        else:
            grd = pystrain.grid.generate_grid(sta_list_ell, args.x_grid_step, args.y_grid_step, True)
        print('[DEBUG] Grid Information:')
        print('[DEBUG]\tLongtitude : from {} to {} with step {} (deg)'.format(grd.x_min, grd.x_max, grd.x_step))
        print('[DEBUG]\tLatitude   : from {} to {} with step {} (deg)'.format(grd.y_min, grd.y_max, grd.y_step))
        print('[DEBUG] Number of Strain Tensors to be estimated: {}'.format(grd.xpts*grd.ypts))
        if fstats:
            print('{:^10s} {:^10s} {:^10s} {:^12s} {:^12s} {:^12s}'.format('Longtitude','Latitude','# stations', 'D (optimal)','CutOff dis.', 'Sigma'), file=fstats)
            print('{:^10s} {:^10s} {:^10s} {:^12s} {:^12s} {:^12s}'.format('deg.','deg.','#', 'Km','#', '/'), file=fstats)
        vprint('[DEBUG] Estimating strain tensor for each cell center:')
        ##  Iterate through the grid (on each cell center). Grid returns cell-centre
        ##+ coordinates in lon/lat pairs, in degrees!
        if args.multiproc_mode:
            grd1, grd2, grd3, grd4 = grd.split2four()
            print('--> grid split to four!')
            fout1=open(".out.thread1", "w")
            fout2=open(".out.thread2", "w")
            fout3=open(".out.thread3", "w")
            fout4=open(".out.thread4", "w")
            if fstats:
                fstats1=open(".sta.thread1", "w")
                fstats2=open(".sta.thread2", "w")
                fstats3=open(".sta.thread3", "w")
                fstats4=open(".sta.thread4", "w")
            else:
                fstats1 = fstats2 = fstats3 = fstats4 = None
            print('[DEBUG] Estimating strain tensors in multi-threading mode')
            #print('--> Thread will be given grid : X:{:}/{:}/{:} Y:{:}/{:}/{:}'.format(grd1.x_min, grd1.x_max, grd1.x_step, grd1.y_min, grd1.y_max, grd1.y_step))
            #print('--> Thread will be given grid : X:{:}/{:}/{:} Y:{:}/{:}/{:}'.format(grd2.x_min, grd2.x_max, grd2.x_step, grd2.y_min, grd2.y_max, grd2.y_step))
            #print('--> Thread will be given grid : X:{:}/{:}/{:} Y:{:}/{:}/{:}'.format(grd3.x_min, grd3.x_max, grd3.x_step, grd3.y_min, grd3.y_max, grd3.y_step))
            #print('--> Thread will be given grid : X:{:}/{:}/{:} Y:{:}/{:}/{:}'.format(grd4.x_min, grd4.x_max, grd4.x_step, grd4.y_min, grd4.y_max, grd4.y_step))
            p1 = multiprocessing.Process(target=compute__, args=(grd1, sta_list_utm, lcm, fout1, fstats1, vprint), kwargs=dargs)
            p2 = multiprocessing.Process(target=compute__, args=(grd2, sta_list_utm, lcm, fout2, fstats2, vprint), kwargs=dargs)
            p3 = multiprocessing.Process(target=compute__, args=(grd3, sta_list_utm, lcm, fout3, fstats3, vprint), kwargs=dargs)
            p4 = multiprocessing.Process(target=compute__, args=(grd4, sta_list_utm, lcm, fout4, fstats4, vprint), kwargs=dargs)
            [ p.start() for p in [p1, p2, p3, p4]]
            [ p.join()  for p in [p1, p2, p3, p4]]
            for fl in [fout1, fout2, fout3, fout4]:
                if not fl.closed:
                    fl.close()
            if fstats:
                for fl in [fstats1, fstats2, fstats3, fstats4]:
                    if not fl.closed:
                        fl.close()
            ##  Note that fout? and fstats? are now closed! We need to
            ##+ concatenate the files though.
            with open(STRAIN_OUT_FILE, 'a') as fout:
                for fnr in range(1,5):
                    with open(".out.thread"+str(fnr), "r") as slave_f:
                        fout.write(slave_f.read())
                    os.remove(".out.thread"+str(fnr))
            if fstats:
                with open(STATISTICS_FILE, 'a') as fstats:
                    for fnr in range(1,5):
                        with open(".sta.thread"+str(fnr), "r") as slave_f:
                            fstats.write(slave_f.read())
                        os.remove(".sta.thread"+str(fnr))
           
        else:
            compute__(grd, sta_list_utm, lcm, fout, fstats, vprint, **dargs)
    else:
        ##  Using veis method. Compute delaunay triangles and estimate one tensor
        ##+ per triangle centre
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
            sstr = ShenStrain(cx, cy, cy<0e0, [sta_list_utm[trng[0]], sta_list_utm[trng[1]], sta_list_utm[trng[2]]], weighting_function='equal_weights')
            sstr.estimate()
            sstr.print_details(fout, lcm)
            ## Print the triangle in the corresponding file (ellipsoidal crd, degrees)
            print('> {:}, {:}, {:}'.format(sta_list_utm[trng[0]].name, sta_list_utm[trng[1]].name, sta_list_utm[trng[2]].name), file=dlnout)
            print('{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}\n{:+8.5f} {:8.5f}'.format(*[ degrees(x) for x in [sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat, sta_list_ell[trng[1]].lon, sta_list_ell[trng[1]].lat, sta_list_ell[trng[2]].lon, sta_list_ell[trng[2]].lat, sta_list_ell[trng[0]].lon, sta_list_ell[trng[0]].lat]]), file=dlnout)
        dlnout.close()
        fout.close()

    ##  Before exiting, write the station information to a file
    write_station_info(sta_list_ell)
    print('[DEBUG] Total running time: {:10.2f} sec.'.format((time.time() - start_time)))
