#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys
import numpy
import operator
from math import atan2, exp, sqrt, floor, pi, degrees
from station import Station
##############################################  pystrain
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *
import pystrain.grid


def barycenter(sta_list):
    ''' Compute the barycenter from a list of stations. The function will use
        each station's self.x and self.y components.
        Barycenter's coordinates will have the same units as the input ones.
    '''
    y_mean = sta_list[0].lat
    x_mean = sta_list[0].lon
    for i in range(1, len(sta_list)):
        y_mean = (sta_list[i].lat + (i-1)*y_mean) / float(i)
        x_mean = (sta_list[i].lon + (i-1)*x_mean) / float(i)
    return y_mean, x_mean

def z_weights(sta_lst, cx, cy, debug_mode=False):
    """ Given a list of Stations and the coordinates of a central point (i.e.
        cx, cy), compute and return the function:
        Z(i) = n*θ(i) / 4π
        which is used as a weighting function fom strain estimation in Shen 
        et al, 2015, see Equation (5a).
        The individual station weights are returned in a list, in the order they
        were passed in in the sta_lst list.
        It is assumed, that the coordinates of every point in the list and the
        cx,cy coordinates, are given in meters.
        Note that, Easting is considered as 'x' or longtitude.
        Northing is considered as 'y' or latitude.

        Parameters:
        -----------
        sta_lst: list
                A list of Station instances. For each one a weight will be
                computed and returned.
        cx:      float
                The x-component of the central point
        cy:      float
                The y-component of the central point

        Returns:
        --------
        list
            Each element in the list is the weight of the respective station in
            the input station list.
    """
    n         = len(sta_lst)
    azimouths = []
    thetas    = []
    #  Get the azimouth of each line from central point (cx,cy) to each point in
    #+ the list. The computed azimouths are stored in a sorted list
    #+ (aka azimouths). This list is made up of dictionary elements, where each
    #+ dictionary contains 1. the azimouth value (in radians) as 'az' and 2. the
    #+ index of the station in the sta_lst, as 'nr'.
    for idx, sta in enumerate(sta_lst):
        az = atan2(sta.lon-cx, sta.lat-cy)
        azimouths.append({'az': az+int(az<0)*2*pi, 'nr': idx}) # normalize to [0, 2pi]
    azimouths = sorted(azimouths, key=operator.itemgetter('az'))
    #  Confirm that all azimouths are in the range [0,2*pi)
    for a in azimouths: assert a['az'] >= 0e0 and a['az'] < 2*pi
    if debug_mode:
        print('[DEBUG] Here are the azimouths:')
        for o in azimouths:
            print('\t{:4s} -> az = {:+8.4f}'.format(sta_lst[o['nr']].name, degrees(o['az'])))
    #  Make a list of the 'theta' angles; for each point, the theta angle is an
    #+ azimouth difference, of the previous minus the next point.
    #  Special care for the first and last elements (theta angles)
    #thetas.append({'w': azimouths[n-1]['az'] - azimouths[1]['az'], 'nr':azimouths[0]['nr']})
    thetas.append({'w': 2e0*pi+(azimouths[1]['az'] - azimouths[n-1]['az']), 'nr':azimouths[0]['nr']})
    for j in range(1, n-1):
        thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
    thetas.append({'w': 2e0*pi+(azimouths[0]['az'] - azimouths[n-2]['az']), 'nr':azimouths[n-1]['nr']})
    if debug_mode:
        print('[DEBUG] Here are the theta angles:')
        for t in thetas:
            print('\t{:4s} -> theta = {:+8.4f}'.format(sta_lst[t['nr']].name, degrees(t['w'])))
    #  Double-check !! All theta angles must be in the range [0, 2*π)
    for angle in thetas:
        assert angle['w'] >= 0 and angle['w'] <= 2*pi, '[ERROR] Error computing statial weights. Station is \"{}\".'.format(sta_lst[angle['nr']].name)
    #  Now compute z = n*θ/4π and re-arrange so that we match the input sta_lst
    if debug_mode:
        print('[DEBUG] Here are the final weights:')
        tmp_lst = [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]
        for idx, val in enumerate(tmp_lst):
            print('\t{:4s} -> z = {:+8.4f}'.format(sta_lst[idx].name, val))
    return [ x['w']*float(n)/(4e0*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]

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
        (as function input) the Z weights.

        Note that the coordinates (both in sta_lst and in cx, cy), must be given
        in meters. Also, cx is matched to sta_lst[i].lon and cy is matched to
        sta_lst[i].lat.

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
                dstep: Step for searching for optimal D, from dmin to dmax until
                    we hit Wt; default is 2
                d:     Value D for computing the weights; if given, then
                       it is used and no effort is made to find an 'optimal'
                       D value.
                debug_mode: Sets debuging mode on.

        Returns:
            tuple: (list, float)
            list is a list of weights (i.e. L(i) values) for each station, in the
            order they are passed in.
            float: is the D value used to compute the weights.

        Raises:
            RuntimeError if dmin > dmax, or if we cannot find an optimal D.
    """
    if 'ltype' not in kargs : kargs['ltype'] = 'gaussian'
    if 'Wt'    not in kargs : kargs['Wt']    = 6
    if 'dmin'  not in kargs : kargs['dmin']  = 1 	
    if 'dmax'  not in kargs : kargs['dmax']  = 500
    if 'dstep' not in kargs : kargs['dstep'] = 1

    debug_mode = False
    if 'debug_mode' in kargs: debug_mode = kargs['debug_mode']

    #  Note: d and dri must be in the same units (here km).
    def gaussian(dri, d):  return exp(-pow(dri/d,2))
    def quadratic(dri, d): return 1e0/(1e0+pow(dri/d,2))

    if kargs['dmin'] >= kargs['dmax'] or kargs['dstep'] < 0:
        raise RuntimeError

    if kargs['ltype'] == 'gaussian':
        l_i = gaussian
    elif kargs['ltype'] == 'quadratic':
        l_i = quadratic
    else:
        raise RuntimeError

    #  Distances for each point from center in km.
    dr = [ sqrt((x.lon-cx)*(x.lon-cx)+(x.lat-cy)*(x.lat-cy))/1000e0 for x in sta_lst ]

    #  If 'd' is given (at input), just compute and return the weights
    if 'd' in kargs:
        d = float(kargs['d'])
        print('[DEBUG] Using passed in \'d\' coef = {:}'.format(d))
        if debug_mode:
            print('[DEBUG] Here are the l weights: (D={:8.5f})'.format(d))
            for i,s in enumerate(sta_lst):
                print('\t{:} Distance: {:5.2f}km., L = {:8.5f}'.format(s.name, dr[i], l_i(dr[i],d)))
        return [ l_i(dri,d) for dri in dr ], d
    #  Else, iterate through [dmin, dmax] to find optimal d; then compute and
    #+ return the weights.
    for d in numpy.arange(kargs['dmin'], kargs['dmax'], kargs['dstep']):
        l    = [ l_i(dri,d) for dri in dr ]
        w    = sum([ x[0]*x[1] for x in zip(l,z_weights) ])*2 # w(i) = l(i)*z(i)
        if int(round(w)) >= kargs['Wt']:
            print('[DEBUG] Found optimal \'d\' coef = {:}'.format(d))
            if debug_mode:
                print('[DEBUG] Here are the l weights: (D={:})'.format(d))
                for i,s in enumerate(sta_lst):
                    print('\t{:} Distance: {:5.2f}km., L = {:8.5f}'.format(s.name, dr[i], l[i]))
            return l, d
    # Fuck! cannot find optimal D
    print('[ERROR] Cannot compute optimal D in weighting scheme')
    raise RuntimeError

def ls_matrices_veis4(sta_lst, cx, cy):
    """  4-parameter deformation
         Dx = dx      + y rotation + x scale
         Dy =      dy - x rotation + y scale
    """
    # numper of rows (observations)
    N = len(sta_lst)*2
    # number of columns (parameters)
    M = 4
    cc  = Station(lon=cx, lat=cy)
    xyr = [ x.distance_from(cc) for x in sta_lst ]
    ## design matrix A, observation matrix b
    A = numpy.zeros(shape=(N,M))
    b = numpy.zeros(shape=(N,1))
    i = 0
    for idx,sta in enumerate(sta_lst):
        dx, dy, dr = xyr[idx]
        A[i]   = [ 1e0, 0e0,  dx, dy ]
        A[i+1] = [ 0e0, 1e0, -dx, dy ]
        b[i]   = (sta.ve/1000) * Wx
        b[i+1] = (sta.vn/1000) * Wy
        i+=2
    assert i is N, "[DEBUG] Failed to construct ls matrices"
    # we can solve this as:
    # numpy.linalg.lstsq(A,b)
    return A, b

def ls_matrices_veis6(sta_lst, cx, cy):
    """  6-parameter deformation
    """
    # numper of rows (observations)
    N = len(sta_lst)*2
    # number of columns (parameters)
    M = 6
    cc  = Station(lon=cx, lat=cy)
    xyr = [ x.distance_from(cc) for x in sta_lst ]
    ## design matrix A, observation matrix b
    A = numpy.zeros(shape=(N,M))
    b = numpy.zeros(shape=(N,1))
    i = 0
    for idx,sta in enumerate(sta_lst):
        dx, dy, dr = xyr[idx]
        A[i]   = [ 1e0, 0e0,  dy,  dx, dy,  0e0]
        A[i+1] = [ 0e0, 1e0, -dx, 0e0, dx,   dy]
        #A[i]    = [ 1e0, 0e0,  dx,  dy, 0e0, 0e0]
        #A[i+1]  = [ 0e0, 1e0, 0e0, 0e0, -dx, dy ]
        b[i]   = sta.ve/1000
        b[i+1] = sta.vn/1000
        i+=2
    assert i is N, "[DEBUG] Failed to construct ls matrices"
    # we can solve this as:
    # numpy.linalg.lstsq(A,b)
    return A, b

def ls_matrices_shen(sta_lst, cx, cy, **kargs):
    """ Construct Least Squares Matrices (A and b) to be solved for. The function
        will first evaluate the covariance matrix C, where:
        W(i) = C(i) * G(i)**(-1), where G(i) = L(i) * Z(i) and C(i) the 1/std.dev
        of each observable.
        Note that we have two observables (x and y) for each station.
    """
    # minimum weight for an observation to be included in the matrices
    min_l_weight = 1e-6
    # numper of rows (observations)
    N = len(sta_lst)*2
    # number of columns (parameters)
    M = 6
    # the (spatial) weights, i.e. Z(i)
    zw = z_weights(sta_lst, cx, cy)
    # the distance weights
    lw, d_coef = l_weights(sta_lst, cx, cy, zw, **kargs)
    assert len(zw) == N/2 and len(zw) == len(lw), '[ERROR] Invalid weight arrays size.'
    # we are only going to use the observations above the minimum L weight threshold
    Nl = sum(1 for l in lw if l >= min_l_weight) * 2
    assert Nl <= N
    # the weight matrix W = Q^(-1) = C^(-1)*G = C^(-1)*(L*Z), which is actualy
    # a column vector
    # W = numpy.zeros(shape=(N,1))
    # distances, dx and dy for each station from (cx, cy). Each element of the
    # array is [ ... (dx, dy, dr) ... ]
    cc  = Station(lon=cx, lat=cy)
    xyr = [ x.distance_from(cc) for x in sta_lst ]
    ## design matrix A, observation matrix b
    A = numpy.zeros(shape=(Nl,M))
    b = numpy.zeros(shape=(Nl,1))
    i = 0
    debug_list = []
    for idx,sta in enumerate(sta_lst):
        sta_used = True
        dx, dy, dr = xyr[idx]
        Wx     = (1e-3/sta.se)*zw[idx]*lw[idx]
        Wy     = (1e-3/sta.sn)*zw[idx]*lw[idx]
        if lw[idx] >= min_l_weight:
            A[i]   = [ Wx*j for j in [1e0, 0e0,  dy,  dx, dy,  0e0] ]
            A[i+1] = [ Wy*j for j in [0e0, 1e0, -dx,   0e0, dx, dy] ]
            b[i]   = sta.ve * Wx
            b[i+1] = sta.vn * Wy
	    i += 2
        else:
            sta_used = False
        debug_list.append({'name':sta.name, 'x':sta.lon, 'y':sta.lat, 'dr':dr/1000e0, 'zw':zw[idx], 'lw':lw[idx], 'wx':Wx, 'wy':Wy, 'used':sta_used})
    assert i == Nl, "[DEBUG] Failed to construct ls matrices"
    # we can solve this as:
    # numpy.linalg.lstsq(A,b)
    for db in sorted(debug_list, key=lambda k: k['dr']):
        print('{:s} {:10.3f}m {:10.3f}m {:10.3f}km z={:10.5f} l={:10.5f} wx={:10.5f} wy={:10.5f} used for strain: {:}'.format(db['name'], db['x'], db['y'], db['dr'],db['zw'],db['lw'],db['wx'],db['wy'],db['used']))
    return A, b

def __strain_info__(str_params):
    info_dict = {}
    tx  = str_params['Ux']
    ty  = str_params['Uy']
    exx = str_params['taux']
    exy = str_params['tauxy'] + str_params['omega']
    eyx = str_params['tauxy'] - str_params['omega']
    eyy = str_params['tauy']
    info_dict['e'] = (exy - eyx) / 2e0
    info_dict['k'] = (exx + eyy) * 1000000. / 2e0
    info_dict['strain'] = sqrt((exx-eyy) * (exx-eyy) + (exy+eyx) * (exy+eyx)) * 1000000.
    info_dict['k_max'] = info_dict['k'] + info_dict['strain'] / 2e0
    info_dict['k_min'] = info_dict['k'] - info_dict['strain'] / 2e0
    info_dict['az'] = degrees(2e0 * atan2(exy + eyx, eyy - exx))
    return info_dict

def __cmp_strain__(str_params, str_params_cov=None):
    info_dict = {}
    x1  = str_params['taux']
    x2  = str_params['tauxy']
    x3  = str_params['tauy']
    cov = pi / 180.0e0
    ##  estimate principle strain rates emax, emin, maximum shear tau_max, 
    ##+ and dextral tau_max azimuth
    emean = (x1+x3) / 2.0e0
    ediff = (x1-x3) / 2.0e0
    taumax= sqrt(x2**2 + ediff**2)
    emax  = emean+taumax
    emin  = emean-taumax
    azim  = -atan2(x2, ediff) / cov / 2.0e0
    #azim  = 90.0e0+azim
    dexazim = azim+45.0e0-180.0e0
    dilat = x1+x3
    return emean, ediff, taumax, emax, emin, azim, dilat

"""
def __write_strain__(strain, fout=None):
    print('longitude latitude      vx+dvx       vy+dvy    cxy         w+dw         exx+dexx         exy+dexy         eyy+deyy         emax+demax       emin+demin        shr+dshr         azi+dazi        dilat+ddilat     dis.    weight    chisq   nsite')
    print('  deg       deg          mm/a         mm/a      /      nano-rad/a     nano-strain/a    nano-strain/a    nano-strain/a    nano-strain/a    nano-strain/a    nano-strain/a        degree       nano-strain/a      km        /         /       /')

print(''.format(strain.x, strain.y, strain.Ux, 0e0, strain.Uy, 0e0, 0e0,
     strain.tauy, 0e0, strain.omega, 0e0, strain.taux, 0e0, strain.tauxy, 0e0,
     emax, 0e0, emin, 0e0, taumax, 0e0, dexazim, 0e0, dilat, 0e0)
"""
class ShenStrain:
    def __init__(self, x=0e0, y=0e0, station_list=[]):
        self.__stalst__ = station_list
        self.__xcmp__   = x
        self.__ycmp__   = y
        self.__zweights__ = None
        self.__lweights__ = None
        self.__options__  = {'ltype': 'gaussian', 'Wt': 24, 'dmin': 1, 'dmax': 500, 'dstep': 2}
        self.__parameters__ = {'Ux':0e0, 'Uy':0e0, 'omega':0e0, 'taux':0e0, 'tauxy':0e0, 'tauy':0e0}
    
    def azimouths(self):
        n         = len(self.__stalst__)
        azimouths = []
        #  Get the azimouth of each line from central point (cx,cy) to each point in
        #+ the list. The computed azimouths are stored in a sorted list
        #+ (aka azimouths). This list is made up of dictionary elements, where each
        #+ dictionary contains 1. the azimouth value (in radians) as 'az' and 2. the
        #+ index of the station in the sta_lst, as 'nr'.
        for idx, sta in enumerate(self.__stalst__):
            az = atan2(sta.lon-self.__xcmp__, sta.lat-self.__ycmp__)
            azimouths.append({'az': az+int(az<0)*2*pi, 'nr': idx}) # normalize to [0, 2pi]
        azimouths = sorted(azimouths, key=operator.itemgetter('az'))
        #  Confirm that all azimouths are in the range [0,2*pi)
        for a in azimouths: assert a['az'] >= 0e0 and a['az'] < 2*pi
        return azimouths

    def azimouth_coverage(self):
        # compute the azimouth to each station; this is a sorted list! azimouth
        # values are in radians.
        azs = self.azimouths()
        n   = len(azs)
        max_angle = degrees(abs(azs[0]['az'] - azs[n-1]['az']))
        return max_angle

    def info(self):
        return __strain_info__(self.__parameters__)

    def print_details(self, fout):
	utm_zone = 34 #added to convert utm 2 latlon
	clat, clon = utm2ell(self.__xcmp__,  self.__ycmp__ , utm_zone) #added to conv utm 2 latlon
        emean, ediff, taumax, emax, emin, azim, dilat =  __cmp_strain__(self.__parameters__)
        print('{:9.5f} {:9.5f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f} {:+10.1f}'.format(degrees(clat), degrees(clon), self.value_of('Ux')*1e3, self.value_of('Uy')*1e3, self.value_of('tauy')*1e9, self.value_of('omega')*1e9, self.value_of('taux')*1e9, self.value_of('tauxy')*1e9, emax*1e9, emin*1e9, taumax*1e9, azim, dilat*1e9), file=fout)

    def set_options(self, **kargs):
        for opt in kargs:
            if opt not in self.__options__:
                print('[DEBUG] Option {:} not relevant for ShenStrain. Ommiting')
            else:
                self.__options__[opt] = kargs[opt]

    def value_of(self, key):
        if key == 'x':
            return self.__xcmp__
        if key == 'y':
            return self.__ycmp__
        if key in self.__parameters__:
            return self.__parameters__[key]
        if key in self.__options__:
            return self.__options__[key]
        raise RuntimeError
    
    def set_xy(self, x, y):
        self.__xcmp__ = x
        self.__ycmp__ = y

    def set_to_barycenter(self):
        self.__ycmp__, self.__xcmp__ = barycenter(self.__stalst__)

    def compute_z_weights(self):
        self.__zweights__ = z_weights(self.__stalst__, self.__xcmp__, self.__ycmp__, debug_mode=False)

    def compute_l_weights(self, **kargs):
        self.__lweights__, d_coef = l_weights(self.__stalst__, self.__xcmp__, self.__ycmp__, self.__zweights__, **self.__parameters__)

    def estimate(self, **kargs):
        """ TODO this is wrong! ls_matrices_shen will re-compute all weights!
        """
        if not self.__zweights__:
            self.compute_z_weights()
        if not self.__lweights__:
            self.compute_l_weights()
        assert len(self.__zweights__) == len(self.__lweights__)
        print('[DEBUG] Info on strain estimation; point at {:10.3f}, {:10.3f}'.format(self.__xcmp__, self.__ycmp__))
        A, b = ls_matrices_shen(self.__stalst__, self.__xcmp__, self.__ycmp__, **kargs)
        m, n = A.shape
        if m <= 3:
            raise RuntimeError('Too few obs to perform LS.')
        #tasp = Station(lon=self.__xcmp__, lat=self.__ycmp__)
        #for idx, sta in enumerate(self.__stalst__):
        #    dx, dy, dr = tasp.distance_from(sta)
        #    print('{:s} {:10.3f}m {:10.3f}m {:10.3f}km z={:} l={:}'.format(sta.name, sta.lon, sta.lat, dr/1000.0, self.__zweights__[idx], self.__lweights__[idx]))
        estim, res, rank, sing_vals = numpy.linalg.lstsq(A, b)
        self.__parameters__['Ux']    = float(estim[0])
        self.__parameters__['Uy']    = float(estim[1])
        self.__parameters__['omega'] = float(estim[2])
        self.__parameters__['taux']  = float(estim[3])
        self.__parameters__['tauxy'] = float(estim[4])
        self.__parameters__['tauy']  = float(estim[5])
        return estim

class VeisStrain:
    def __init__(self, x=0e0, y=0e0, station_list=[]):
        self.__stalst__ = station_list
        self.__xcmp__   = x
        self.__ycmp__   = y
        self.__options__  = {}
        self.__parameters__ = {'Ux':0e0, 'Uy':0e0, 'omega':0e0, 'taux':0e0, 'tauxy':0e0, 'tauy':0e0}
        # self.__parameters__ = {'tx':0e0, 'ty':0e0, 'exx':0e0, 'exy':0e0, 'eyx':0e0, 'eyy':0e0}

    def value_of(self, key):
        if key == 'x':
            return self.__xcmp__
        if key == 'y':
            return self.__ycmp__
        if key in self.__parameters__:
            return self.__parameters__[key]
        raise RuntimeError
    
    def info(self):
        return __strain_info__(self.__parameters__)

    def set_xy(self, x, y):
        self.__xcmp__ = x
        self.__ycmp__ = y

    def set_to_barycenter(self):
        self.__ycmp__, self.__xcmp__ = barycenter(self.__stalst__)

    def estimate(self,parameters=6):
        if parameters == 4:
            A, b = ls_matrices_veis4(self.__stalst__, self.__xcmp__, self.__ycmp__)
        elif parameters == 6:
            A, b = ls_matrices_veis6(self.__stalst__, self.__xcmp__, self.__ycmp__)
        else:
            raise RuntimeError
        estim, res, rank, sing_vals = numpy.linalg.lstsq(A, b)
        if parameters == 6:
            self.__parameters__['Ux']    = float(estim[0])
            self.__parameters__['Uy']    = float(estim[1])
            self.__parameters__['omega'] = float(estim[2])
            self.__parameters__['taux']  = float(estim[3])
            self.__parameters__['tauxy'] = float(estim[4])
            self.__parameters__['tauy']  = float(estim[5])
        else:
            raise RuntimeError
        return estim

    def info(self):
        info_dict = {}
        tx  = self.__parameters__['Ux']
        ty  = self.__parameters__['Uy']
        exx = self.__parameters__['taux']
        exy = self.__parameters__['tauxy'] + self.__parameters__['omega']
        eyx = self.__parameters__['tauxy'] - self.__parameters__['omega']
        eyy = self.__parameters__['tauy']
        info_dict['e'] = (exy - eyx) / 2e0
        info_dict['k'] = (exx + eyy) * 1000000. / 2e0
        info_dict['strain'] = sqrt((exx-eyy) * (exx-eyy) + (exy+eyx) * (exy+eyx)) * 1000000.
        info_dict['k_max'] = info_dict['k'] + info_dict['strain'] / 2e0
        info_dict['k_min'] = info_dict['k'] - info_dict['strain'] / 2e0
        info_dict['az'] = degrees(2e0 * atan2(exy + eyx, eyy - exx))
        return info_dict
