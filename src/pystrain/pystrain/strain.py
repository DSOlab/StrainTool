#! /usr/bin/python

import os, numpy, operator
from math import atan2
from station import Station

class Grid:
    def __init__(self, x_min, x_max, x_step, y_min, y_max, y_step):
        self.y_min = y_min
        self.y_max = y_max
        self.y_step= y_step
        self.x_min = x_min
        self.x_max = x_max
        self.x_step= x_step
        self.cy    = self.y_min + self.y_step/2
        self.cx    = self.x_min + self.x_step/2
        self.iter    = 0

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

def make_grid(sta_list, y_step, x_step):
    y_min = sys.float_info.max
    y_max = sys.float_info.min
    x_min = sys.float_info.max
    x_max = sys.float_info.min
    for s in sta_lst:
        if   s.x > max_x:
            max_x = s.x
        elif s.x < min_x:
            min_x = s.x
        if   s.y > max_y:
            max_y = s.y
        elif s.y < min_y:
            min_y = s.y
    # adjust max and min to step
    div = (y_max-y_min)/y_step
    rem = div - math.floor(div)
    y_min -= rem/2
    y_max += rem/2
    div = (x_max-x_min)/x_step
    rem = div - math.floor(div)
    x_min -= rem/2
    x_max += rem/2
    # return a Grid instance
    return Grid(y_min, y_max, y_step, x_min, x_max, x_step)

def z_weights(sta_lst, cx, cy):
    n = len(sta_lst)
    azimouths = []
    thetas    = []
    i = 0
    for sta in sta_lst:
        az = atan2(sta.lon-cx, sta.lat-cy)
        azimouths.append({'az': az+int(az<0)*2*pi, 'nr': i}) # normalize to [0, 2pi]
        i+=1
    azimouths = sorted(azimouths, key=operator.itemgetter('az'))
    thetas.append({'w': azimouths[1]['az'] - azimouths[n-1]['az'], 'nr':azimouths[0]['nr']})
    for j in range(1, n-1):
        thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
    thetas.append({'w':azimouths[0]['az'] - azimouths[n-2]['az'] + 2*pi, 'nr':azimouths[n-1]['nr']})
    return [ x['w']*n/4*pi for x in sorted(thetas, key=operator.itemgetter('nr') ]

"""
def ls_matrices(sta_lst):
    ## numper of rows (observations)
    N = len(sta_lst)*2
    ## number of columns (parameters)
    M = 6
    ## the barycenter as Station
    blat, blon = barycenter(sta_lst)
    baryc = Station(lat=blat, lon=blon)
    ## design matrix A, observation matrix b
    i = 0
    A = numpy.zeros(shape=(N,M))
    b = numpy.zeros(shape=(M,1))
    for sta in sta_lst:
        dlon, dlat, r = baryc.distance_from(sta)
        A[i]   = [1, 0, dlat,  dlon, dlat, 0   ]
        A[i+1] = [0, 1, -dlon, 0,    dlon, dlat]
        b[i]   = sta.ve
        b[i+1] = sta.vn
        i+=2
"""

if __name__ == "__main__":
    g = Grid(23.123, 39.2432, 3, 20.123, 47.123, 5)
    for f,l in g:
        print f,l
