#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys
import numpy
from scipy import linalg
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

def __cmp_strain__(str_params, params_cov=None):
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
    sec_inv = sqrt(x1*x1+x2*x2+x3*x3)
    #if params_cov != None:
    # cut the part of vcv that holds tau* info
    vcv = params_cov[3:6, 3:6]
    # estimate sigma of tau_max
    v = numpy.zeros(shape=(3,1))
    v[0,:] = (x1-x3)/4.0e0/taumax
    v[1,:] = x2/taumax
    v[2,:] = -v[0,:]
    staumax = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
    # estimate sigma of emax
    v[0,:] = 0.5e0*(1+(x1-x3)/2.e0/taumax)
    v[1,:] = x2/taumax
    v[2,:] = 0.5e0*(1-(x1-x3)/2.e0/taumax)
    semax = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
    # estimate sigma of emin
    v[0,:] = 0.5e0*(1-(x1-x3)/2.0e0/taumax)
    v[1,:] = -x2/taumax
    v[2,:] = 0.5e0*(1+(x1-x3)/2.0e0/taumax)
    semin = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
    # estimate sigma of azimuth
    cf = 1.0e0/((x1-x3)**2e0+4.0e0*x2**2e0)
    v[0,:] = cf*x2
    v[1,:] = -cf*(x1-x3)
    v[2,:] = -v[0,:]
    sazim = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
    # estimate sigma of dilatation
    sdilat = sqrt(vcv[0,0]+vcv[2,2]+2e0*vcv[0,2])
    return emean, ediff, taumax, staumax, emax, semax, emin, semin, azim, sazim, dilat, sdilat, sec_inv

class ShenStrain:
    def __init__(self, x=0e0, y=0e0, station_list=[], **kwargs):
        self.__stalst__ = station_list
        self.__xcmp__   = x
        self.__ycmp__   = y
        self.__zweights__ = None
        self.__lweights__ = None
        self.__options__    = {'ltype': 'gaussian', 'Wt': 24, 'dmin': 1, 'dmax': 500, 'dstep': 2, 'd_coef': None, 'cutoff_dis': None}
        self.__parameters__ = {'Ux':0e0, 'Uy':0e0, 'omega':0e0, 'taux':0e0, 'tauxy':0e0, 'tauy':0e0}
        self.__vcv__ = None
        for key in kwargs:
            if key in self.__options__:
                self.__options__[key] = kwargs[key]
            #else:
            #    print('[DEBUG] Irrelevant key in Strain constructor: {:}; skipped'.format(key))
        if self.__options__['ltype'] == 'gaussian':
            self.__options__['cutoff_dis'] = 2.15e0
        else:
            self.__options__['cutoff_dis'] = 10e0

    def clean_weight_matrices(self):
        self.__zweights__ = None
        self.__lweights__ = None

    def filter_sta_wrt_distance(self, d=None):
        cc = Station(lon=self.__xcmp__, lat=self.__ycmp__)
        if not d: d = self.__options__['d_coef']
        limit = self.__options__['cutoff_dis'] * d
        nlst = [ s for s in self.__stalst__ if s.distance_from(cc)[2] <= limit*1e3 ]
        return nlst
    
    def azimouths(self, other_sta_lst=None):
        #  Get the azimouth of each line from central point (cx,cy) to each point in
        #+ the list. The computed azimouths are stored in a sorted list
        #+ (aka azimouths). This list is made up of dictionary elements, where each
        #+ dictionary contains:
        #+   1. the azimouth value (in radians) as 'az' and
        #+   2. the index of the station in the sta_lst, as 'nr'.
        #  To compute the azimouths, we need the ΔX and ΔY components, which are
        #+ computed as:
        #+ ΔX = sta.lon-self.__xcmp__
        #+ ΔY = sta.lat-self.__ycmp__
        #+ Obviously, all of these quantities **should** be in a cartesian RF.
        #  For example, if the instance's station list is:
        #+ [ankr, buku, dion, ...] and the function returns: 
        #+ [{'az':0.34, 'nr':2}, ..., {'az':3.01, 'nr':0}]
        #+ it means, that the line from the instance's centre to dion has an 
        #+ azimouth of 0.34 radians, the line from the instance's centre to 
        #+ ankr has an azimouth of 3.01 radians, etc...
        #  Note that sorted means **in ascending order** (obviously).
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        azimouths = []
        for idx, sta in enumerate(stalst):
            az = atan2(sta.lon-self.__xcmp__, sta.lat-self.__ycmp__)
            azimouths.append({'az': az+int(az<0)*2*pi, 'nr': idx}) # normalize to [0, 2pi]
        azimouths = sorted(azimouths, key=operator.itemgetter('az'))
        #  Confirm that all azimouths are in the range [0,2*pi)
        for a in azimouths: assert a['az'] >= 0e0 and a['az'] < 2*pi
        return azimouths

    def ls_matrices(self, **kargs):
        """ Construct Least Squares Matrices (A and b) to be solved for. The function
            will first evaluate the covariance matrix C, where:
            W(i) = C(i) * G(i)**(-1), where G(i) = L(i) * Z(i) and C(i) the 1/std.dev
            of each observable.
            Note that we have two observables (x and y) for each station.
        """
        # numper of rows (observations)
        N = len(self.__stalst__)*2
        # number of columns (parameters)
        M = 6
        # the (spatial) weights, i.e. Z(i)
        if not self.__zweights__ or not self.__lweights__:
            raise RuntimeError
        d_coef = self.__options__['d_coef']
        zw     = self.__zweights__
        lw     = self.__lweights__
        assert len(zw) == N/2 and len(zw) == len(lw), '[ERROR] Invalid weight arrays size.'
        # the weight matrix W = Q^(-1) = C^(-1)*G = C^(-1)*(L*Z), which is actualy
        # a column vector
        # W = numpy.zeros(shape=(N,1))
        # distances, dx and dy for each station from (cx, cy). Each element of the
        # array is [ ... (dx, dy, dr) ... ]
        cc  = Station(lon=self.__xcmp__, lat=self.__ycmp__)
        xyr = [ x.distance_from(cc) for x in self.__stalst__ ]
        ## design matrix A, observation matrix b
        A = numpy.zeros(shape=(N,M))
        b = numpy.zeros(shape=(N,1))
        i = 0
        #debug_list = []
        for idx, sta in enumerate(self.__stalst__):
            dx, dy, dr = xyr[idx]
            Wx     = (1e-3/sta.se)*sqrt(zw[idx]*lw[idx])
            Wy     = (1e-3/sta.sn)*sqrt(zw[idx]*lw[idx])
            A[i]   = [ Wx*j for j in [1e0, 0e0,  dy,  dx, dy,  0e0] ]
            A[i+1] = [ Wy*j for j in [0e0, 1e0, -dx,   0e0, dx, dy] ]
            b[i]   = sta.ve * Wx
            b[i+1] = sta.vn * Wy
            i += 2
            #debug_list.append({'name':sta.name, 'x':sta.lon, 'y':sta.lat, 'dr':dr/1000e0, 'zw':zw[idx], 'lw':lw[idx], 'wx':Wx, 'wy':Wy, 'used':sta_used})
        assert i == N, "[DEBUG] Failed to construct ls matrices"
        # we can solve this as:
        # numpy.linalg.lstsq(A,b)
        #for db in sorted(debug_list, key=lambda k: k['dr']):
        #    print('{:s} {:10.3f}m {:10.3f}m {:10.3f}km z={:10.5f} l={:10.5f} wx={:10.5f} wy={:10.5f} used for strain: {:}'.format(db['name'], db['x'], db['y'], db['dr'],db['zw'],db['lw'],db['wx'],db['wy'],db['used']))
        return A, b

    def z_weights(self, other_sta_lst=None, debug_mode=False):
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
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        thetas = self.compute_theta_angles(stalst)
        assert len(thetas) == n
        #  Now compute z = n*θ/4π and re-arrange so that we match the input sta_lst
        if debug_mode:
            print('[DEBUG] Here are the final weights:')
            tmp_lst = [ x['w']*n/(4*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]
            for idx, val in enumerate(tmp_lst):
                print('\t{:4s} -> z = {:+8.4f}'.format(stalst[idx].name, val))
        # return [ x['w']*float(n)/(4e0*pi) for x in sorted(thetas, key=operator.itemgetter('nr')) ]
        wt_az = 0.25e0
        azi_avrg = wt_az * 360.0e0 / n
        azi_tot = (1.0e0+wt_az)*360.0e0
        # wght(i)=(0.5d0*(daz1-daz2)+azi_avrg) *nslct/azi_tot ! compute azimuthal weighting
        return [ (0.5e0*degrees(x['w'])+azi_avrg)*n/azi_tot for x in sorted(thetas, key=operator.itemgetter('nr')) ]

    def compute_theta_angles(self, other_sta_lst=None):
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        azimouths = self.azimouths(stalst)
        assert len(azimouths) == n
        thetas = []
        #  Make a list of the 'theta' angles; for each point i, the theta angle
        #+ is an azimouth difference, of the previous (i-1) minus the next 
        #+ point (i+1).
        #  Special care for the first and last elements (theta angles).
        thetas.append({'w': 2e0*pi+(azimouths[1]['az'] - azimouths[n-1]['az']), 'nr':azimouths[0]['nr']})
        for j in range(1, n-1):
            thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'], 'nr':azimouths[j]['nr']})
        thetas.append({'w': 2e0*pi+(azimouths[0]['az'] - azimouths[n-2]['az']), 'nr':azimouths[n-1]['nr']})
        #  Double-check !! All theta angles must be in the range [0, 2*π)
        for angle in thetas:
            assert angle['w'] >= 0 and angle['w'] <= 2*pi, '[ERROR] Error computing statial weights. Station is \"{}\".'.format(stalst[angle['nr']].name)
        return thetas

    def l_weights(self, other_sta_lst=None, **kargs):
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
        debug_mode = False
        if 'debug_mode' in kargs: debug_mode = kargs['debug_mode']

        #  Note: d and dri must be in the same units (here km).
        def gaussian(dri, d):  return exp(-pow(dri/d,2))
        def quadratic(dri, d): return 1e0/(1e0+pow(dri/d,2))

        if self.__options__['ltype'] == 'gaussian':
            l_i = gaussian
        elif self.__options__['ltype'] == 'quadratic':
            l_i = quadratic
        else:
            raise RuntimeError
        
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst

        #  Distances for each point from center in km.
        dr = [ sqrt((x.lon-self.__xcmp__)*(x.lon-self.__xcmp__)+(x.lat-self.__ycmp__)*(x.lat-self.__ycmp__))/1000e0 for x in stalst ]
        d  = float(self.__options__['d_coef'])
        return [ l_i(dri,d) for dri in dr ], d

    def find_optimal_d(self):
        assert self.__options__['dmin'] < self.__options__['dmax']
        for d in numpy.arange(self.__options__['dmin'], self.__options__['dmax'], self.__options__['dstep']):
            self.__options__['d_coef'] = d
            new_sta_lst = self.filter_sta_wrt_distance(d)
            if len(new_sta_lst) > 3:
                lwghts,_ = self.l_weights(new_sta_lst)
                zwghts   = self.z_weights(new_sta_lst)
                assert len(lwghts) == len(zwghts)
                w = sum([ x[0]*x[1] for x in zip(lwghts,zwghts) ])*2 # w(i) = l(i)*z(i)
                if int(round(w)) >= self.__options__['Wt']:
                    print('[DEBUG] Found optimal D value {}; Number of stations involved {} out of {}.'.format(d, len(new_sta_lst), len(self.__stalst__)))
                    return lwghts, zwghts, d
        # Fuck! cannot find optimal D
        print('[ERROR] Cannot compute optimal D in weighting scheme')
        raise RuntimeError

    def beta_angles(self):
        azimouths = self.azimouths()
        n         = len(azimouths)
        betas     = []
        #  Make a list of the 'beta' angles; for each point, the theta angle is an
        #+ azimouth difference, of the previous minus the next point.
        #  Special care for the first and last elements (beta angles)
        betas.append(2e0*pi+(azimouths[0]['az'] - azimouths[n-1]['az']))
        for j in range(0, n-2):
            betas.append(azimouths[j+1]['az'] - azimouths[j]['az'])
        betas.append(2e0*pi+(azimouths[0]['az'] - azimouths[n-1]['az']))
        #  Double-check !! All theta angles must be in the range [0, 2*π)
        for angle in betas:
            assert angle >= 0 and angle <= 2*pi
        assert len(betas) == n
        return betas

    def info(self):
        return __strain_info__(self.__parameters__)

    def print_details(self, fout):
	utm_zone = 34 #added to convert utm 2 latlon
	clat, clon = utm2ell(self.__xcmp__,  self.__ycmp__ , utm_zone) #added to conv utm 2 latlon
        emean, ediff, taumax, staumax, emax, semax, emin, semin, azim, sazim, dilat, sdilat, sec_inv =  __cmp_strain__(self.__parameters__, self.__vcv__)
        print('{:9.5f} {:9.5f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f} {:+7.1f}'.format(degrees(clat), degrees(clon), self.value_of('Ux')*1e3, sqrt(self.__vcv__[0,0])*1e3, self.value_of('Uy')*1e3, sqrt(self.__vcv__[1,1])*1e3, self.value_of('omega')*1e9, sqrt(self.__vcv__[2,2])*1e9, self.value_of('taux')*1e9, sqrt(self.__vcv__[3,3])*1e9, self.value_of('tauxy')*1e9, sqrt(self.__vcv__[4,4])*1e9, self.value_of('tauy')*1e9, sqrt(self.__vcv__[5,5])*1e9, emax*1e9, semax*1e9, emin*1e9, semin*1e9, taumax*1e9, staumax*1e9, azim, sazim, dilat*1e9, sdilat*1e9, sec_inv*1e9), file=fout)

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

    def estimate(self, **kargs):
        if not self.__options__['d_coef']:
            print('[DEBUG] Searching for optimal D parameter.')
            if self.__options__['dmin'] >= self.__options__['dmax'] or self.__options__['dstep'] < 0:
                raise RuntimeError
            lwghts, zwghts, d = self.find_optimal_d()
            self.__stalst__ = self.filter_sta_wrt_distance(d)
        else:
            d = self.__options__['d_coef']
            print('[DEBUG] Using optimal D parameter {}km.'.format(d))
            self.__stalst__ = self.filter_sta_wrt_distance(d)
            lwghts,_ = self.l_weights()
            zwghts   = self.z_weights()
        self.__zweights__ = zwghts
        self.__lweights__ = lwghts
        A, b = self.ls_matrices(**kargs)
        VcV  = numpy.dot(A.T, A)
        m, n = A.shape
        if m <= 3:
            raise RuntimeError('Too few obs to perform LS.')
        ##  Note: To silence warning in versions > 1.14.0, use a third argument,
        ##+ rcond=None; see https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
        estim, res, rank, sing_vals = numpy.linalg.lstsq(A, b)
        # Parameter variance-covariance matrix
        try:
            # A-posteriori variance
            sigma0_post = float(res[0])
            print('[DEBUG] A-posteriori std. deviation = {:}'.format(sqrt(sigma0_post)))
            bvar = linalg.inv(VcV) * sigma0_post
            self.__vcv__ = bvar
        except:
            print('[DEBUG] Cannot compute var-covar matrix! Probably singular.')
            self.__vcv__ = None
        self.__parameters__['Ux']    = float(estim[0])
        self.__parameters__['Uy']    = float(estim[1])
        self.__parameters__['omega'] = float(estim[2])
        self.__parameters__['taux']  = float(estim[3])
        self.__parameters__['tauxy'] = float(estim[4])
        self.__parameters__['tauy']  = float(estim[5])
        if self.__vcv__ is None:
            raise ArithmeticError('Failed to Compute var-covar matrix of parameters')
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
