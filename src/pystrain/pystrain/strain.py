#! /usr/bin/python

import sys, numpy, operator
from math import atan2, exp, sqrt, floor, pi, degrees
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
        self.x_min = x_min
        self.x_max = x_max
        self.x_step= x_step
        self.y_min = y_min
        self.y_max = y_max
        self.y_step= y_step
        self.cxi    = 0
        self.cyi    = 0
        self.xpts   = (x_max-x_min)/x_step - 1
        self.ypts   = (y_max-y_min)/y_step - 1

    def __iter__(self):
        return self

    def xidx2xval(self, idx):
        return self.x_min + self.x_step/2 + self.x_step*idx
    
    def yidx2yval(self, idx):
        return self.y_min + self.y_step/2 + self.y_step*idx

    def next(self):
        if self.cxi > self.xpts:
            if self.cyi > self.ypts:
                raise StopIteration
            self.cxi  = 0
            self.cyi += 1
            return self.x_max - self.x_step/2, self.yidx2yval(self.cyi-1)
        else:
            self.cxi += 1
            return self.xidx2xval(self.cxi-1), self.yidx2yval(self.cyi)

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
    print "\tRegion: Easting: {}/{} Northing: {}/{}".format(x_min, x_max, y_min, y_max)
    s      = float((floor((y_max-y_min)/y_step)+1)*y_step)
    r      = s-(y_max-y_min)
    y_min -= r/2
    y_max += r/2
    print "\tAdjusted Northing: from {} to {} with step={} pts={}".format(y_min, y_max, y_step, (y_max-y_min)/y_step)
    s      = float((floor((x_max-x_min)/x_step)+1)*x_step)
    r      = s-(x_max-x_min)
    x_min -= r/2
    x_max += r/2
    print "\tAdjusted Easting: from {} to {} with step={} pts={}".format(x_min, x_max, x_step, (x_max-x_min)/x_step)
    # return a Grid instance
    return Grid(x_min, x_max, x_step, y_min, y_max, y_step)

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
    #print '\t[DEBUG] Station azimouths (decimal degrees):'
    #for ii in azimouths:
    #    print '\t\tstation {} A={}'.format(sta_lst[ii['nr']].name, degrees(ii['az']))
    thetas.append({'w': azimouths[1]['az'] - azimouths[n-1]['az'] + 2*pi, 'nr':azimouths[0]['nr']})
    for j in range(1, n-1):
        thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
    thetas.append({'w':azimouths[0]['az'] - azimouths[n-2]['az'] + 2*pi, 'nr':azimouths[n-1]['nr']})
    # double check!
    for angle in thetas:
        assert angle['w'] >= 0 and angle['w'] <= 2*pi, '[ERROR] Error computing statial weights. Station is \"{}\".'.format(sta_lst[angle['nr']].name)
    #print '\t[DEBUG] Station theta angles (decimal degrees):'
    #for ii in thetas:
    #    print '\t\tstation {} theta={}'.format(sta_lst[ii['nr']].name, degrees(ii['w']))
    #return [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]
    return [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]
    #for z in sorted(thetas, key=operator.itemgetter('nr')):
    #    #print '\t[DEBUG] {} Z={}'.format(sta_lst[z['nr']].name, z['w']*n/(4*pi))
    #return Z

def l_weights(sta_lst, cx, cy, z_weights, Wt=24, dmin=1, dmax=500, dstep=2):
    """
        dmin, dmax and dstep must be given in km.
    """
    if dmin >= dmax or dstep < 0:
        raise RuntimeError
    # distance from center for each station
    dr = [ sqrt((x.lon-cx)*(x.lon-cx)+(x.lat-cy)*(x.lat-cy))/1000 for x in sta_lst ]
    # iterate through [dmin, dmax] to find optimal d
    for d in numpy.arange(dmin, dmax, dstep):
        #print '\t[DEBUG] Computing Li for d={}'.format(d)
        l    = [ exp(-pow(dri/d,2)) for dri in dr ]
        w    = sum([ x[0]*x[1] for x in zip(l,z_weights) ])*2 # w(i) = l(i)*z(i)
        #for idx,li in enumerate(l):
        #    print '\t\t{} L={} (distance={}km, value={})'.format(sta_lst[idx].name, li, dr[idx], l[idx])
        #print '\t[DEBUG] Sum of weights = {}'.format(w)
        if int(round(w)) == Wt:
            return l
        #print '\t[DEBUG] l_weights: w={} != Wt={}'.format(w, Wt)
    # fuck! cannot find optimal D
    print '[ERROR] Cannot compute optimal D in weighting scheme'
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
