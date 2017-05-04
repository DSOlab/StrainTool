#! /usr/bin/python

import sys, numpy, operator
from math import atan2, exp, sqrt, floor, pi
from station import Station

class Grid:
    """
        A very simple Grid class to be used within the StrainTensor project.
        A Grid instance has x- and y- axis limits and step sizes (i.e x_min,
        x_max, x_step, y_min, y_max, y_step).
        It is iterable; when iterating, the instance will return the center
        of each cell, starting from the bottom left corner and ending at the
        top right. The iteration is performed row-wise (i.e.
        > [x0+xstep/2, y0+ystep/2]
        > [x0+xstep/2, y0+ystep/2+ystep]
        > [x0+xstep/2, y0+ystep/2+2*ystep]
        > ...
        > [x0+xstep/2, ymax-ystep/2]
        > [x0+xstep/2+xstep, y0+ystep/2]
        > [x0+xstep/2+xstep, y0+ystep/2+ystep]
        > [x0+xstep/2+xstep, y0+ystep/2+2*ystep]

    """
    def __init__(self, x_min, x_max, x_step, y_min, y_max, y_step):
        self.y_min = y_min
        self.y_max = y_max
        self.y_step= y_step
        self.x_min = x_min
        self.x_max = x_max
        self.x_step= x_step
        self.cy    = self.y_min + self.y_step/2
        self.cx    = self.x_min + self.x_step/2

    def __iter__(self):
        return self

    def next(self):
        if self.cx >=  self.x_max - self.x_step/2:
            if self.cy >= self.y_max - self.y_step/2:
                raise StopIteration
            self.cx =  self.x_min + self.x_step/2
            self.cy += self.y_step
            return self.cy - self.y_step, self.x_max - self.x_step/2
        else:
            self.cx += self.x_step
            return self.cy, self.cx - self.x_step

def barycenter(sta_list):
    '''
        Compute the barycenter from a list of stations
    '''
    y_mean = sta_list[0].y
    x_mean = sta_list[0].x
    for i in range(1, len(sta_list)):
        y_mean = (sta_list[i].y + (i-1)*y_mean)/i
        x_mean = (sta_list[i].x + (i-1)*x_mean)/i
    return y_mean, x_mean

def make_grid(sta_lst, x_step, y_step):
    """
        Given a list of Stations and x- and y-axis step sizes, compute and
        return a Grid instance. The max/min x and y values (of the Grid) are
        extracted from the station coordinates; if needed, they are adjusted
        so that (xmax-xmin) is divisible (without remainder) with xstep.
    """
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
    div = (y_max-y_min)/y_step
    rem = div - floor(div)
    y_min -= rem/2
    y_max += rem/2
    div = (x_max-x_min)/x_step
    rem = div - floor(div)
    x_min -= rem/2
    x_max += rem/2
    # return a Grid instance
    return Grid(y_min, y_max, y_step, x_min, x_max, x_step)

def z_weights(sta_lst, cx, cy):
    """
        Given a  list of Stations and the coordinates of a central point (i.e.
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
    thetas.append({'w': azimouths[1]['az'] - azimouths[n-1]['az'], 'nr':azimouths[0]['nr']})
    for j in range(1, n-1):
        thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
    thetas.append({'w':azimouths[0]['az'] - azimouths[n-2]['az'] + 2*pi, 'nr':azimouths[n-1]['nr']})
    return [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]

def l_weights(sta_lst, cx, cy, z_weights, Wt=24, dmin=1, dmax=100, dstep=1):
    """
    """
    if dmin >= dmax or dstep < 0:
        raise RuntimeError
    # distance from center for each station
    dr = [ sqrt((x.lon-cx)*(x.lon-cx)+(x.lat-cy)*(x.lat-cy)) for x in sta_lst ]
    # iterate through [dmin, dmax] to find optimal d
    for d in numpy.arange(dmin, dmax, dstep):
        l    = [ exp(-pow(dri/d,2)) for dri in dr ]
        w    = sum([ x[0]*x[1] for x in zip(l,z_weights) ]) # w(i) = l(i)*z(i)
        if w == Wt: return l
        print '\t[DEBUG] l_weights: w={} != Wt={}'.format(w, Wt)
    # fuck! cannot find optimal D
    print '[DEBUG] Cannot compute optimal D in weighting scheme'
    raise RuntimeError

def ls_matrices(sta_lst, cx, cy):
    # numper of rows (observations)
    N = len(sta_lst)*2
    # number of columns (parameters)
    M = 6
    # the (spatial) weights, i.e. Z(i)
    zw = z_weights(sta_lst, cx, cy)
    # the distance weights
    lw = l_weights(sta_lst, cx, cy, zw)
    # the weight matrix W = Q^(-1) = C^(-1)*G = C^(-1)*(L*Z), which is actualy
    # a column vector
    # W = numpy.zeros(shape=(N,1))
    # distances, dx and dy for each station from (cx, cy). Each element of the
    # array is [ ... (dx, dy, dr) ... ]
    cc  = Station(lon=cx, lat=cy)
    xyr = [ x.distance_from(cc) for x in sta_lst ]
    ## design matrix A, observation matrix b
    A = numpy.zeros(shape=(N,M))
    b = numpy.zeros(shape=(M,1))
    i = 0
    for sta in sta_lst:
        dx, dy, dr = xyr[i]
        Wx     = (1.0/sta.ve)*zw[i]*lw[i]
        Wy     = (1.0/sta.vn)*zw[i]*lw[i]
        A[i]   = [1, 0,  dy,  dx, dy,  0] * Wx
        A[i+1] = [0, 1, -dx,   0, dx, dy] * Wy
        b[i]   = sta.ve * Wx
        b[i+1] = sta.vn * Wy
        i+=2
    assert i is N*2, "[DEBUG] Failed to construct ls matrices"
    # we can solve this as:
    # numpy.linalg.lstsq(A,b)
    return A, b

if __name__ == "__main__":
    g = Grid(23.123, 39.2432, 3, 20.123, 47.123, 5)
    for f,l in g:
        print f,l
