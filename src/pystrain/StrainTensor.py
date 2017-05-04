#! /usr/bin/python

## standard libls
import sys
from copy import deepcopy
from math import degrees, floor
## numpy
import numpy
## pystrain
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *

X_GRID_STEP = 40000
Y_GRID_STEP = 40000
print '[DEBUG] x and y step sizes are {} and {} meters on UTM'.format(X_GRID_STEP, Y_GRID_STEP)

##  Parse stations from input file
sta_list_ell = parse_ascii_input( sys.argv[1] )
print '\t[DEBUG] Number of stations parsed: {}'.format(len(sta_list_ell))

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

##  Construct the grid, based on station coordinates (Ref. UTM)
grd = make_grid(sta_list_utm, X_GRID_STEP, Y_GRID_STEP)
print '[DEBUG] Constructed the grid. Limits are:'
print '\t[DEBUG] Easting : from {} to {} with step {}'.format(grd.x_min, grd.x_max, grd.x_step)
print '\t[DEBUG] Northing: from {} to {} with step {}'.format(grd.y_min, grd.y_max, grd.y_step)

print '[DEBUG] Estimating strain tensor for each cell center'
##  Iterate through the grid (on each cell center)
for x, y in grd:
    print '\t[DEBUG] Strain at E={}, N={}'.format(x, y)
    ##  Construct the LS matrices, A and b
    A, b = ls_matrices(sta_list_utm, x, y)
    ##  Solve LS
    x, res, rank, sing_vals = numpy.linalg.lstsq(A, b)
    ## print result
    print '\tUx={}\n\tUy={}\n\tomega={}\n\tTx={}\n\tTxy={}\n\tTy={}'.format(x[0], x[1], x[2], x[3], x[4], x[5])
