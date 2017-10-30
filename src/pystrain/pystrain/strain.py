#! /usr/bin/python
#-*- coding: utf-8 -*-

import sys, numpy, operator
from math import atan2, exp, sqrt, floor, pi, degrees

from station import Station
from grid    import Grid

def barycenter(sta_list):
    ''' Compute the barycenter from a list of stations
    '''
    y_mean = sta_list[0].y
    x_mean = sta_list[0].x
    for i in range(1, len(sta_list)):
        y_mean = (sta_list[i].y + (i-1)*y_mean) / float(i)
        x_mean = (sta_list[i].x + (i-1)*x_mean) / float(i)
    return y_mean, x_mean

def make_grid(sta_lst, x_step, y_step):
    """ Given a list of Stations and x- and y-axis step sizes, compute and
        return a Grid instance. The max/min x and y values (of the Grid) are
        extracted from the station coordinates; if needed, they are adjusted
        so that (xmax-xmin) is divisible (without remainder) with xstep.
    """
    print "[DEBUG] Constructing grid in UTM"
    y_min = sys.float_info.max
    y_max = sys.float_info.min
    x_min = sys.float_info.max
    x_max = sys.float_info.min
    for s in sta_lst:
        if   s.lon > x_max:
            x_max = s.lon
        elif s.lon < x_min:
            x_min = s.lon
        if   s.lat > y_max:
            y_max = s.lat
        elif s.lat < y_min:
           y_min = s.lat
    # adjust max and min to step
    print "\t[DEBUG] Region: Easting: {}/{} Northing: {}/{}".format(x_min, x_max, y_min, y_max)
    s      = float((floor((y_max-y_min)/y_step)+1)*y_step)
    r      = s-(y_max-y_min)
    y_min -= r/2
    y_max += r/2
    print "\t[DEBUG] Adjusted Northing: from {} to {} with step={} pts={}".format(y_min, y_max, y_step, (y_max-y_min)/y_step)
    s      = float((floor((x_max-x_min)/x_step)+1)*x_step)
    r      = s-(x_max-x_min)
    x_min -= r/2
    x_max += r/2
    print "\t[DEBUG] Adjusted Easting: from {} to {} with step={} pts={}".format(x_min, x_max, x_step, (x_max-x_min)/x_step)
    # return a Grid instance
    return Grid(x_min, x_max, x_step, y_min, y_max, y_step)

def z_weights(sta_lst, cx, cy):
    """ Given a list of Stations and the coordinates of a central point (i.e.
        cx, cy), compute and return the function:
        Z(i) = n*theta(i)/4pi
        which is used as a weighting function from strain estimation in Shen 
        et al, 2015, see Equation (5a).
        The individual station weights are returned in a list, in the order they
        were passed in in the sta_lst list.
    """
    n         = len(sta_lst)
    azimouths = []
    thetas    = []
    for idx, sta in enumerate(sta_lst):
        az = atan2(sta.lon-cx, sta.lat-cy)
        azimouths.append({'az': az+int(az<0)*2*pi, 'nr': idx}) # normalize to [0, 2pi]
    azimouths = sorted(azimouths, key=operator.itemgetter('az'))
    thetas.append({'w': azimouths[1]['az'] - azimouths[n-1]['az'] + 2*pi, 'nr':azimouths[0]['nr']})
    for j in range(1, n-1):
        thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
    thetas.append({'w':azimouths[0]['az'] - azimouths[n-2]['az'] + 2*pi, 'nr':azimouths[n-1]['nr']})
    # double check!
    for angle in thetas:
        assert angle['w'] >= 0 and angle['w'] <= 2*pi, '[ERROR] Error computing statial weights. Station is \"{}\".'.format(sta_lst[angle['nr']].name)
    return [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]

def l_weights(sta_lst, cx, cy, z_weights, **kargs):
    """ Compute L(i) for each of the points in the station list sta_lst, where
        L(i) = exp(-ΔR(i)**2/D**2) -- Gaussian, or
        L(i) = 1/(1+ΔR(i)**2/D**2) -- Quadratic

        To compute the function, the parameter D is needed; in this function the
        approach is to try different D values (from dmin to dmax with step
        dstep), untill the total re-weighting coefficient of the data approaches
        the limit Wt. The re-weighting coefficient is:
        W = Σ{i=0, i=len(sta_lst)*2}G(i)
        where Σ denotes the sum and G(i) = L(i) * Z(i); this is why we also need
        (as function input) the Z weights (as list).

        The parameters dmin, dmax and dstep must be given in km. Wt must be an
        integer.

        Args:
            sta_lst: A list of Station instances, i.e. the list of stations
            cx: The x coordinate of the center point
            cy: The y coordinate of the center point
            z_weights: A list of weights (i.e. float values) for each station
            **kwargs: Arbitrary keyword arguments; the following are relevant:
                ltype: type of weighting function; can be either 'gaussian', or
                    'quadratic'.
                Wt:    value for optimal Wt; default is 24
                dmin:  Minimum value of tested D's; default is 1
                dmax:  Maximum value of tested D's; defult is 500
                dstep: Step for searching for optimal D, from dmin to dmax untill
                    we hit Wt; default is 2

        Returns:
            a list of weights (i.e. L(i) values) for each station.

        Raises:
            RuntimeError if dmin > dmax, or if we cannot find an optimal D.
    """
    if not 'ltype' in kargs : kargs['ltype'] = 'gaussian'
    if not 'Wt' in kargs    : kargs['Wt']    = 24
    if not 'dmin' in kargs  : kargs['dmin']  = 1
    if not 'dmax' in kargs  : kargs['dmax']  = 500
    if not 'dstep' in kargs : kargs['dstep'] = 2

    def gaussian(dri, d):  return exp(-pow(dri/d,2))
    def quadratic(dri, d): return 1.0/(1.0+pow(dri/d,2))

    if kargs['dmin'] >= kargs['dmax'] or kargs['dstep'] < 0: raise RuntimeError

    if kargs['ltype'] == 'gaussian':
        l_i = gaussian
    elif kargs['ltype'] == 'quadratic':
        l_i = quadratic
    else:
        raise RuntimeError

    # list of distances for each point from center
    dr = [ sqrt((x.lon-cx)*(x.lon-cx)+(x.lat-cy)*(x.lat-cy))/1000 for x in sta_lst ]
    # iterate through [dmin, dmax] to find optimal d
    for d in numpy.arange(kargs['dmin'], kargs['dmax'], kargs['dstep']):
        #l    = [ exp(-pow(dri/d,2)) for dri in dr ]
        l    = [ l_i(dri,d) for dri in dr ]
        w    = sum([ x[0]*x[1] for x in zip(l,z_weights) ])*2 # w(i) = l(i)*z(i)
        if int(round(w)) == kargs['Wt']:
            return l
    # fuck! cannot find optimal D
    print '[ERROR] Cannot compute optimal D in weighting scheme'
    raise RuntimeError

def ls_matrices(sta_lst, cx, cy, **kargs):
    """ Construct Least Squares Matrices (A and b) to be solved for. The function
        will first evaluate the covariance matrix C, where:
        W(i) = C(i) * G(i)**(-1), where G(i) = L(i) * Z(i) and C(i) the 1/std.dev
        of each obervable.
    """
    # numper of rows (observations)
    N = len(sta_lst)*2
    # number of columns (parameters)
    M = 6
    # the (spatial) weights, i.e. Z(i)
    zw = z_weights(sta_lst, cx, cy)
    # the distance weights
    lw = l_weights(sta_lst, cx, cy, zw, **kargs)
    assert len(zw) == N/2 and len(zw) == len(lw), '[ERROR] Invalid weight arrays size.'
    # the weight matrix W = Q^(-1) = C^(-1)*G = C^(-1)*(L*Z), which is actualy
    # a column vector
    # W = numpy.zeros(shape=(N,1))
    # distances, dx and dy for each station from (cx, cy). Each element of the
    # array is [ ... (dx, dy, dr) ... ]
    cc  = Station(lon=cx, lat=cy)
    xyr = [ x.distance_from(cc) for x in sta_lst ]
    ## design matrix A, observation matrix b
    A = numpy.zeros(shape=(N,M))
    b = numpy.zeros(shape=(N,1))
    i = 0
    for idx,sta in enumerate(sta_lst):
        dx, dy, dr = xyr[idx]
        Wx     = (1.0/sta.ve)*zw[idx]*lw[idx]
        Wy     = (1.0/sta.vn)*zw[idx]*lw[idx]
        A[i]   = [ Wx*j for j in [1, 0,  dy,  dx, dy,  0] ]
        A[i+1] = [ Wy*j for j in [0, 1, -dx,   0, dx, dy] ]
        b[i]   = (sta.ve/1000) * Wx
        b[i+1] = (sta.vn/1000) * Wy
        i+=2
    assert i is N, "[DEBUG] Failed to construct ls matrices"
    # we can solve this as:
    # numpy.linalg.lstsq(A,b)
    return A, b

if __name__ == "__main__":
    g = Grid(23.123, 39.2432, 3, 20.123, 47.123, 5)
    for f,l in g:
        print f,l
