#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys
import numpy
from scipy import linalg
import operator
from math import atan2, exp, sqrt, floor, pi, degrees

import pystrain.grid
from pystrain.station import Station
from pystrain.strain import *
from pystrain.geodesy.utm import *
from pystrain.iotools.iparser import *

DEBUG_MODE = False

def barycenter(sta_list):
    ''' Compute the barycenter from a list of stations. The function will use
        each station's self.lat and self.lon components.
        Barycenter's coordinates will have the same units as the input ones.

        Args:
            sta_lst (list of Station): a list of Station instances.

        Returns:
            tuple (float, float): first element is the average of the stations
                                  longtitude, while the second is the average
                                  latitude.
    '''
    if len(sta_list) == 0:
        raise ValueError("[ERROR] Cannot compute barycentre for empty list of stations")
    y_mean = sta_list[0].lat
    x_mean = sta_list[0].lon
    for i in range(1, len(sta_list)):
        y_mean = (sta_list[i].lat + (i-1)*y_mean) / float(i)
        x_mean = (sta_list[i].lon + (i-1)*x_mean) / float(i)
    return x_mean, y_mean

class ShenStrain:
    """A class to represeent Strain Tensors.

        ToDo: Need to check the following documentation!!

        Attributes:
            __xcmp__ (float): x coordinate (could also be easting)
            __ycmp__ (float): y coordinate (could also be northing)
            __stalst__ (list of Station): a list of stations to be used
                for strain estimation. Some of them may be filtered out.
            __zweights__ (list): a list of floats, where each element is the
                weight of the station based on azimouthal coverage.
            __lweights__ (list): a list of floats, where each element is the
                weight of the station based on distance.
            __options__ (dictionary): A dictionary holding the following:
                * ltype (str): gaussian or quadratic; this is the function
                               to be used for distance weight computation.
                * Wt (float): parameter for finding optimal D value (km)
                * dmin (float): min D value for finding optimal D (km)
                * dmax (float): max D value for finding optimal D (km)
                * dstep (float): step for range dmin, dmax (km)
                * d_coef (float): optimal D value (km)
                * cutoff_dis (float): cut off distance (km)
                * weighting_function (str): can be shen (to use shen weighting
                    algorithm), or equal_weights (to use equal weights for all
                    stations)
                * verbose_mode (bool): sets verbose mde on if True; i.e. print
                  debugging messages
            vprint (function): if the instance is created with with verbose_mode
                on, then this function is just print(); else vprint is a noop.

        Note:
            The Station class has no x or y members; it only has lon and lat.
            Because a lot of (member) functions will need to compute distances
            between points and the Strain Tensor centre point (aka __xcmp__ and
            __ycmp__), it will be assumed that:
                * the Station.lon member corresponds to __xcmp__, and
                * the Station.lat member corresponds to __ycmp__
            For example, to compute the distance between the point with index i
            and the Strain Tensor centre, the instance will do something like:
            d = sqrt((__stalst__[i].lon - __xcmp__)^2 + (__stalst__[i].lat - __ycmp__)^2)
    """

    def __init__(self, x=0e0, y=0e0, in_southern_hemisphere=False, station_list=[], **kwargs):
        """ ShenStrain constructor.

            Args:
                x (float): x coordinate (could also be easting)
                y (float): y coordinate (could also be northing)
                station_list (list of Station): a list of stations to be used
                           for strain estimation. Some of them may be filtered out.
                **kwargs: a dictionary containing any of the keys:
                    * ltype (str): 'gaussian' or 'quadratic'; this is the
                                   function to be used for distance weight
                                   computation.
                    * Wt (float): parameter for finding optimal D value (km)
                    * dmin (float): min D value for finding optimal D (km)
                    * dmax (float): max D value for finding optimal D (km)
                    * dstep (float): step for range dmin, dmax (km)
                    * d_coef (float): optimal D value (km)
                    * cutoff_dis (float): cut off distance (km) --see Warning--
                    * weighting_function (str): can be 'shen' (to use shen 
                      weighting algorithm), or 'equal_weights' (to use equal 
                      weights for all stations)
                    * verbose_mode (bool): sets verbose mode on if True; i.e.
                      print debugging messages.

            Warning:
                the value of __options__[cutoff_dis] will be automatically set.
                Even if the user specifies a value at the constrcutor, it will
                change.
        """
        ##  Set input values (x, y, station_list) and initiallize all others to
        ##+ default values.
        self.__stalst__ = station_list
        self.__xcmp__   = x
        self.__ycmp__   = y
        self.__in_shemisphere__ = in_southern_hemisphere
        self.__zweights__ = None
        self.__lweights__ = None
        self.__options__  = {
            'ltype': 'gaussian',
            'Wt': 24, 
            'dmin': 1,
            'dmax': 500,
            'dstep': 2,
            'd_coef': None,
            'cutoff_dis': None,
            'weighting_function': 'shen',
            'verbose_mode': False
        }
        self.__parameters__ = {
            'Ux':0e0,
            'Uy':0e0,
            'omega':0e0, 
            'taux':0e0,
            'tauxy':0e0,
            'tauy':0e0
        }
        self.__vcv__ = None
        ## Resolve the dictionary passed in (if any)
        for key in kwargs:
            if key in self.__options__:
                self.__options__[key] = kwargs[key]
        if self.__options__['ltype'] == 'gaussian':
            self.__options__['cutoff_dis'] = 2.15e0
        else:
            self.__options__['cutoff_dis'] = 10e0
        ##  If in verbose_mode, set the vprint function to print; else vprint
        ##+ is a noop
        self.vprint = print if self.__options__['verbose_mode'] else lambda *a, **k: None

    def clean_weight_matrices(self):
        """ Set both distance and spatial weight lists (aka zweights and 
            lweigts) to None.
        """
        self.__zweights__ = None
        self.__lweights__ = None

    def filter_sta_wrt_distance(self, d=None):
        """ Filter instance's station list wrt the cutof distance.
            
            This function will compute a cut-off distance:
                COD = self.__options__['cutoff_dis'] * D
            and filter all stations in the instance's __stalst__ list. A new
            station list will be returned, containing only stations that are
            less (or equal to) COD km apart from the instance's centre (aka
            __xcmp__, __ycmp__) point.
            If the input parameter d is provided, then the function will use
            that as D, else the function will use self.__options__['d_coef']
            as D (in COD computation).

            Args:
                d (float): parameter for computing the cut-off distance. If not
                           provided, __options__['d_coef'] will be used. d 
                           should be provided in km.

            Returns:
                list of Station: a list of Station instances, where each station
                                 is less than cut-off distance away from the
                                 instance.
            Note:
                The returned list is not assigned to the instance's __stalst__.
                If you want that, then do it manually.

        """
        if not self.__options__['cutoff_dis']:
            raise ValueError("[ERROR] Cannot filter station list; cutoff_dis is None!")
        cc = Station(lon=self.__xcmp__, lat=self.__ycmp__)
        if not d: d = self.__options__['d_coef']
        limit = self.__options__['cutoff_dis'] * d
        ##  OPT try optimized squared distance (aka remove the square roots).
        ##+ That is instead of filtering based on sqrt(Δx^2 + Δy^2) < limit*1e3
        ##+ we will use (Δx^2 + Δy^2) < limit*limit
        nlst = [ s for s in self.__stalst__ if s.squared_distance_from(cc) <= limit*limit ]
        ## In debug mode, check that we have the correct results
        if DEBUG_MODE:
            nlst1 = [ s for s in self.__stalst__ if s.distance_from(cc)[2] <= limit*1e3 ]
            assert len(nlst) == len(nlst1)
            for i,s in enumerate(nlst1):
                assert s.name == nlst[i].name
        return nlst
    
    def azimouths(self, other_sta_lst=None):
        """ Azimouth of line containing the instance's centre and each point.
            
            Get the azimouth of each line from central point (__xcmp__, __ycmp__) 
            to each point in the other_sta_lst list. The computed azimouths 
            are stored in a sorted list (aka azimouths). This list is made up 
            of dictionary elements, where each dictionary contains:
              1. the azimouth value (in radians) as 'az' and
              2. the index of the station in the other_sta_lst, as 'nr'.
            To compute the azimouths, we need the ΔX and ΔY components, which 
            are computed as:
                ΔX = sta.lon-self.__xcmp__
                ΔY = sta.lat-self.__ycmp__
                az = atan2(ΔX, ΔY), normalized to [0, 2π)
            Obviously, all of these quantities **should** be in a cartesian RF.
            For example, if the instance's station list is:
            [ankr, buku, dion, ...] and the function returns: 
            [{'az':0.34, 'nr':2}, ..., {'az':3.01, 'nr':0}]
            it means, that the line from the instance's centre to dion has an 
            azimouth of 0.34 radians, the line from the instance's centre to 
            ankr has an azimouth of 3.01 radians, etc...
            Note that sorted means **in ascending order** (obviously).
            If the user does not provide a other_sta_lst parameter, then
            the function will use the isntance's __stalst__.

            Args:
                other_sta_lst (list of Station): a list of Station instances.
                    If not provided, __stalst__ will be used.

           Returns:
                (sorted) list: This list is made up of dictionary elements, 
                where each dictionary contains:
                    1. the azimouth value (in radians) as 'az' and
                    2. the index of the station in the other_sta_lst, as 'nr'.

            Note:
                All azimouths will fall in range [0, 2π)
        """
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        azimouths = []
        for idx, sta in enumerate(stalst):
            az = atan2(sta.lon-self.__xcmp__, sta.lat-self.__ycmp__)
            azimouths.append({'az': az+int(az<0e0)*(2e0*pi), 'nr': idx})
        azimouths = sorted(azimouths, key=operator.itemgetter('az'))
        if DEBUG_MODE:
            ##  if in debug mode, confirm that all azimouths are in the range
            ##+ [0,2*pi)
            for a in azimouths: assert a['az'] >= 0e0 and a['az'] < 2*pi
        return azimouths

    def ls_matrices(self, sigma0=1):
        """ Construct Least Squares Matrices (A and b) to be solved for.
        
            Matrix A is the "design matrix" and b is the observation vector.
            This function, will first compute the weight matrix (W). The
            formulation of the weight matrix depends on the option 
            "weighting_function" of the instance (shen or equal_weights).
            Note that the computation of the weight matrix is actually performed
            via a call to this->make_weight_matrix().

            Note that the function make_weight_matrix does not actually return
            the weight matrix but its square root, aka W = P^(1/2). This
            is the matrix we will use here (we only need the square root of
            the weight matrix to multiply with matrices A and b). Also note
            that the weight matrix is not scaled by any a-priori sigma0; hence
            the scaling should be done here (if needed).

            Then the function will form two matrices, A and b. The matrix A,
            is the design matrix multiplied by W, that is each row of A, is:
            A(i,:)   = (1, 0, Δx, Δy, 0, Δy)*W(i)
            A(i+1,:) = (0, 1, 0, Δx, Δy, -Δx)*W(i+1),
            where:
            Δx is the distance (x-component) between station i and the 
               instance's centre in meters
            Δy is the distance (y-component) between station i and the 
               instance's centre in meters
            W(i) is the (scalar) square root of weight of station i, for the x 
               (or east) component
            W(i+1) is the (scalar) square root of weight of station i, for the 
               y (or north) component
            For the matrix (actually vector) b, the function will form:
            b(i)   = sta.ve * W(i)
            b(i+1) = sta.vn * W(i+1)
            
            If we solve A\b, we' ll end up with the parameter vector:
            [ Ux, Uy, τx, τxy, τy, ω ]**T
            The number of rows in A and b matrices, will be:
            len(self.__stalst__)*2
            The units of the paramter vector are:
            [ m, m, strain/yr, strain/yr, strain/yr, strain/yr ]**T

            Args:
                sigma0 (float): A-priori sigma0 (σ0) for the formulation of the
                    weight matrix.

            Returns:
                tuple (numpy.array, numpy.array): first matrix is A, second is
                    b

        """
        ## number of rows (observations)
        N = len(self.__stalst__)*2
        ## number of columns (parameters)
        M = 6
        ## the weights, i.e. σ0 * W(i)
        W = sigma0 * self.make_weight_matrix()
        assert W.shape == (N,1)
        ##  Distances, dx and dy for each station from (cx, cy). Each element
        ##+ of the array is xyr = [ ... (dx, dy, dr) ... ]
        cc  = Station(lon=self.__xcmp__, lat=self.__ycmp__)
        xyr = [ cc.distance_from(x) for x in self.__stalst__ ]
        ## design matrix A, observation matrix b
        A = numpy.zeros(shape=(N,M))
        b = numpy.zeros(shape=(N,1))
        i = 0
        for idx, sta in enumerate(self.__stalst__):
            dx, dy, dr = xyr[idx]
            Wx     = W[i]
            Wy     = W[i+1]
            A[i]   = [ Wx*j for j in [1e0, 0e0,  dx, dy, 0e0, dy] ]
            A[i+1] = [ Wy*j for j in [0e0, 1e0, 0e0, dx, dy, -dx] ]
            b[i]   = sta.ve * Wx
            b[i+1] = sta.vn * Wy
            i += 2
        assert i == N, "[DEBUG] Failed to construct ls matrices"
        return A, b

    def make_weight_matrix(self):
        """ Construct the square root of weight matrix W <- P^(1/2)

            This function will construct the weight matrix to be used for
            strain estimation (via LSE). The weight matrix will be formed
            according to the option in __options['weighting_function']__.
            Note that the returned matrix is actually an array(!!) of shape
            (N,1), where N = len(self.__stalst__)*2. Each element in the
            returned matrix, will be the weight for the corresponding station's
            lon and lat component. E.g., if __stalst__ = ['dyng', 'ankr', ...],
            then W[0] is the x-weight of dyng, W[1] is the y-weight of dyng,
            W[2] is the x-weight of ankr, W[3] is the y-weight of ankr ....

            Note that no a-priori sigma0 (σ0) scaling factor is used here. If
            you want, you can scale the resulting matrix later on.

            If weighting scheme is 'shen', then the function will:
                * read z- and l- weights from the instance's __zweights__ and
                  __lweights__ lists
                * compute for each station a pair of sigmas, as:
                  w_x = (1/σe) * sqrt(Z(i)*L(i))
                  w_y = (1/σn) * sqrt(Z(i)*L(i))
            If weighting scheme is 'equal weights', then the function will
            return a weight matrix with all elements equal to 1.

            Returns:
                numpy matrix (Nx1): The weight matrix computed, of size 
                (2*len(__stalst__),1), where
                W[0] <- weight of __stalst__[0] x component
                W[1] <- weight of __stalst__[0] y component
                W[2] <- weight of __stalst__[1] x component
                [ ... ]
            
            Note:
                If using the 'shen' weighting scheme, then the L and Z weights
                should have already been computed.
        """
        ## number of rows (observations)
        N = len(self.__stalst__)*2
        W = numpy.ones(shape=(N,1))
        ## Use Shen's weighting scheme
        if self.__options__['weighting_function'] == 'shen':
            if not self.__zweights__ or not self.__lweights__:
                raise RuntimeError("[ERROR] Z or L weights not set; cannot compute weight matrices")
            d_coef = self.__options__['d_coef']
            zw = self.__zweights__
            lw = self.__lweights__
            i = 0
            for idx, sta in enumerate(self.__stalst__):
                W[i]   = (1e0/sta.se)*sqrt(zw[idx]*lw[idx])
                W[i+1] = (1e0/sta.sn)*sqrt(zw[idx]*lw[idx])
                i += 2
            assert i == N
        elif self.__options__['weighting_function'] == 'equal_weights':
            self.vprint('[DEBUG] Using equal-weight covar matrix!')
            #pass
        else:
            raise RuntimeError("[ERROR] Invalid weighting function option")
        return W

    def z_weights(self, other_sta_lst=None):
        """ Compute spatial (i.e. azimouthal coverage) weights, accordin to
            Shen et al, 2015
        
            Given a list of Stations, compute and return the function:
            Z(i) = n*θ(i) / 4π
            which is used as a weighting function fom strain estimation in Shen 
            et al, 2015, see Equation (5a).
            The individual station weights are returned in a list, in the
            order they were passed in in the other_sta_lst list.
            It is assumed, that the coordinates of every point in the list and
            the __xcmp__, __ycmp__coordinates, are given in a Cartesian RF and
            are in meters.
            Note that (for each station), Station.lon is considered the 'x'
            component and Station.lat is considered the 'y' component.

            Args:
                other_sta_lst (list of Station) : A list of Station instances. For
                    each one a weight will be computed and returned. If not given,
                    the instance's __stalst__ list will be used.

            Returns:
                list of floats: Each element in the list is the weight of the
                    respective station in the input station list (aka 
                    len(list) == len(other_sta_lst))

            Warning:
                The weighting function is NOT Z(i) = n*θ(i) / 4π, but it is
                extracted from VISR. So actualy, the formula is:
                Z(i) = (.5 * θ(deg.) + azi_avrg) * n / azi_tot,
                where azi_avrg = 0.25 * 360 / n
                and azi_tot = (1+0.25)*3600
        """
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        thetas = self.compute_theta_angles(stalst)
        assert len(thetas) == n
        wt_az = 0.25e0
        azi_avrg = wt_az * 360e0 / n
        azi_tot = (1e0+wt_az)*360e0
        return [ (0.5e0*degrees(a)+azi_avrg)*n/azi_tot for a in thetas ]

    def compute_theta_angles(self, other_sta_lst=None):
        """ Compute θ angles, aka next minus the previous point.

            Make a list of the 'theta' angles; for each point i, the theta angle
            is an azimouth difference, of the previous (i-1) minus the next 
            point (i+1). All of the θ angles will be in the range (0, 2π).
               
               P(i-1)
                \
                 \           P(i)
                  \         /
                   \  β1   /
                    \     /
                     \   /
                      \ /
                       C  β2                    θ = β1 + β2
                        \
                         \
                          \
                           \
                           P(i+1)
            
            Args:
                other_sta_lst: A list of Station instances. For each one a θ (theta)
                    angle is computed and returned. If not given, the instance's
                    __stalst__ list will be used.

            Returns:
                list of floats. Each element in the list, is the θ angle of
                    the corresponding station in other_sta_lst.
        """
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst
        n = len(stalst)
        azimouths = self.azimouths(stalst)
        assert len(azimouths) == n
        thetas = []
        ##  Special care for the first and last elements (theta angles).
        thetas.append({'w' : 2e0*pi+(azimouths[1]['az'] - azimouths[n-1]['az']),\
                       'nr': azimouths[0]['nr']})
        for j in range(1, n-1):
            thetas.append({'w':azimouths[j+1]['az'] - azimouths[j-1]['az'],\
                          'nr':azimouths[j]['nr']})
        thetas.append({'w': 2e0*pi+(azimouths[0]['az'] - azimouths[n-2]['az']),\
                      'nr':azimouths[n-1]['nr']})
        ##  Double-check !! All theta angles must be in the range [0, 2*π)
        if DEBUG_MODE:
            for angle in thetas: assert angle['w'] >= 0 and angle['w'] <= 2*pi
        return [ i['w'] for i in sorted(thetas, key=operator.itemgetter('nr')) ]

    def l_weights(self, other_sta_lst=None):
        """ Compute distance-dependent weights.
        
            Compute L(i) for each of the points in the station list sta_lst,
            where
                L(i) = exp(-ΔR(i)**2/D**2) -- Gaussian, or
                L(i) = 1/(1+ΔR(i)**2/D**2) -- Quadratic
            It is assumed, that the coordinates of every point in the list and
            the __xcmp__, __ycmp__coordinates, are given in a Cartesian RF and
            are in meters.
            Note that (for each station), Station.lon is considered the 'x'
            component and Station.lat is considered the 'y' component.
            The function will extract the gaussian or quadratic function,
            depending on __options__['ltype']. Also the (optimal) D value
            will be extracted from __options__['d_coef'], which should be in km.

            Args:
                other_sta_lst: A list of Station instances, i.e. the list of
                    stations. If not given, the instance's __stalst__ will be
                    used.

            Returns:
                tuple: (list, float)
                list is a list of weights (i.e. L(i) values) for each station,
                in the order they are passed in.
                float is the D value used to compute the weights.

        """
        ##  Note: d and dri must be in the same units (here km).
        def gaussian(dri, d):  return exp(-pow(dri/d,2))
        def quadratic(dri, d): return 1e0/(1e0+pow(dri/d,2))
        ## assign the correct weighting formula
        if self.__options__['ltype'] == 'gaussian':
            l_i = gaussian
        elif self.__options__['ltype'] == 'quadratic':
            l_i = quadratic
        else:
            raise RuntimeError("[ERROR] Invalid distance-dependent weighting function")
        
        stalst = self.__stalst__ if other_sta_lst is None else other_sta_lst

        #  Distances for each point from center in km.
        dr = [ sqrt((x.lon-self.__xcmp__)*(x.lon-self.__xcmp__)+\
                    (x.lat-self.__ycmp__)*(x.lat-self.__ycmp__))/1000e0 \
                    for x in stalst ]
        d = float(self.__options__['d_coef'])
        if not d:
            raise RuntimeError("[ERROR] D-coefficient ton set; cannot compute distance-dependent weights")
        return [ l_i(dri,d) for dri in dr ], d

    def find_optimal_d(self):
        """ Find optimal D coefficient, for distance weighting.

            This function will test the range [dmin, dstep) with a step of
            dstep (km) for an optimal D coefficient (aka the coefficient used
            to compute the distance-dependent weighting).
            The optimal D is that for which:
            int(W) >= int(W_t), where:
            W = Σ{L(i)*Z(i)} for i in 1,2,...,len(stations)*2
            Note that when testing the individual D coefficients, the stations
            list used may not be the instance's __stalst__. That is because,
            for every D selected, the __stalst__ is filtered using:
            self.filter_sta_wrt_distance(D).

            In summary, the function will:
                - loop through D=[dmin:dmax:dstep]
                    * make new station list (newstalst) filtering the instance's
                      station list, using D
                    * compute Z-weights
                    * compute L-weights
                    * compute W = Sum{l(i)*z(i)} for i=1, ..., len(newstalst)
                    * if int(W) >= Wt stop
            
            Note:
                The function will use the instance's:
                    __options__['dmin']
                    __options__['dmax']
                    __options__['dstep']
                    __options__['Wt']
                and set the instance's parameter:
                    __options__['dcoef']

            Returns:
                a tuple with elements:
                -- list (floats): the lweights (i.e. distance weights) computed 
                    with the optimal D coeff.
                -- list (floats): the zweights (i.e. spatial weights) computed
                    with the optimal D coeff
                -- float: The optimal D coefficient in km.

            Warning:
                The returned lists (lweights and zweights), may not be of the
                same size as __stalst__. They actually correspond to the list
                of stations that where filtered using the function:
                self.filter_sta_wrt_distance(D). So if you do something like:
                lw, zw, D = self.find_optimal_d()
                valid_sta = self.filter_sta_wrt_distance(D)
                then, lw[i] and zw[i] will be the weights of valid_sta[i] (and
                not self.__stalst__[i]).

            Raises:
                RuntimeError if no  optimal D coeff can be found within the
                range [dmin, dmax).
        """
        assert self.__options__['dmin'] < self.__options__['dmax']
        for d in numpy.arange(self.__options__['dmin'], \
                              self.__options__['dmax'], \
                              self.__options__['dstep']):
            self.__options__['d_coef'] = d
            new_sta_lst = self.filter_sta_wrt_distance(d)
            if len(new_sta_lst) > 3:
                lwghts,_ = self.l_weights(new_sta_lst)
                zwghts   = self.z_weights(new_sta_lst)
                assert len(lwghts) == len(zwghts)
                w = sum([ x[0]*x[1] for x in zip(lwghts,zwghts) ])*2 # w(i) = l(i)*z(i)
                if int(round(w)) >= int(self.__options__['Wt']):
                    return lwghts, zwghts, d
        # Fuck! cannot find optimal D
        self.vprint('[ERROR] Cannot compute optimal D in weighting scheme')
        raise RuntimeError

    def beta_angles(self):
        """ Return the β angles (internal angles).

            Make a list of the β ('beta') angles; for each point, the β angle
            is an azimouth difference, of the next minus the previous point
            aka β = Az(i+1) - Az(i)

               P(i)
                \
                 \           P(i+1)
                  \         /
                   \       /
                    \  β  /
                     \   /
                      \ /
                       C

            Returns:
                list (float): The β angles (all in range [0, 2π)). The length
                    of the returned list, is the length of __stalst__ minus one.
                    Note that the β angles are **not in correspondance** with 
                    the __stalst__ list (aka, β[0] is not the angle between
                    __stalst__[0] and __stalst__[1]).
        """
        ##  Note that this->azimouths() returns a dictionary, sorted on the
        ##+ the azimouth value. Hence, azimouths[0]['az'] may not be the
        ##+ azimouth of the line from the centre to station[0].
        azimouths = self.azimouths()
        n = len(azimouths)
        betas = []
        betas.append(2e0*pi+(azimouths[0]['az'] - azimouths[n-1]['az']))
        for j in range(0, n-1):
            betas.append(azimouths[j+1]['az'] - azimouths[j]['az'])
        ##  Double-check !! All theta angles must be in the range [0, 2*π)
        if DEBUG_MODE:
            for angle in betas: assert angle >= 0 and angle <= 2*pi
        assert len(betas) == n
        return betas

    def cmp_strain(self, params_cov=None):
        ''' Compute strain tensor parameters and sigmas

            This function will compute the strain tensor parameters:
                * emean
                * ediff
                * taumax
                * emax
                * emin
                * azim
                * dilat
                * sec_inv
            given that the "fundamental" parameters have already been estimated
            (i.e. the parameters [Ux, Uy, τx, τxy, τy, ω] with units:
            [ m, m, strain/yr, strain/yr, strain/yr, strain/yr ]**T ).

            The function will also compute the parameter corresponding std.
            deviation  values (aka staumax, semax, semin, sazim, sdilat), if
            the Variance-Covariance matrix of the "fundamental" parameters is
            passed in (as params_cov).

            The variance-covariance matrix, must be a numpy.array of size
            6x6. It is assumed, that the rows concering the τx, τxy, τy 
            parameters, are in indexes [2,5), i.e. the var-covar matrix is of
            type:
            | σ_Ux^2  σ_UxUy  σ_Uxτx σ_Uxτxy  σ_Uxτy   σ_Uxω  | (row 0)
            | σUxUy   σ_Uy^2  σ_Uyτx σ_Uyτxy  σ_Uyτy   σ_Uyω  | (row 1)
            | σUxτx   .       σ_τx^2 σ_τxτxy  σ_τxτy   σ_τxω  | (row 2)
            | .       .       .      σ_τxy^2  σ_τxyτy  σ_τxyω | (row 3)
            | .       .       .      .        σ_τy^2   σ_τyω  | (row 4)
            | .       .       .      .        .        σ_ω^2  | (row 5)
        col:  0       1       2      3        4         5       

            Actualy, to perform the calculation, the function is going to cut
            the submatrix params_cov[2:5, 2:5], so all other elements could
            contain whatever values.
            If the function parameter params_cov is set to None, then all values
            for the parameters std. deviations, will be also set to None.

            Args:
                params_cov (numpy.matrix(6,6)): The variance-covariance matrix
                of the parameters [Ux, Uy, τx, τxy, τy, ω] --in that order--.

            Returns:
                a tuple, holding the elements:
                (emean, ediff, taumax, staumax, emax, semax, emin, semin, \
                azim, sazim, dilat, sdilat, sec_inv, ssec_inv)
                If params_cov is None, then the values of the sigmas (staumax,
                semax, semin, sazim, sdilat and ssec_inv) are all set to 'None'.
                All units are strain/year

            Note:
                Normaly, the user should call the funtion with self.__vcv__ as
                the (input) params_cov parameter.

            Refs:
                The functions to compute the strain parameters, are taken from
                Shen's VISR fortran code.
        '''
        x1  = self.__parameters__['taux']   ## strain/yr
        x2  = self.__parameters__['tauxy']  ## strain/yr
        x3  = self.__parameters__['tauy']   ## strain/yr
        cov = pi / 180e0
        ##  estimate principle strain rates emax, emin, maximum shear tau_max, 
        ##+ and dextral tau_max azimuth
        emean = (x1+x3) / 2e0               ## strain/yr
        ediff = (x1-x3) / 2e0               ## strain/yr
        taumax= sqrt(x2**2 + ediff**2)      ## strain/yr
        emax  = emean+taumax                ## strain/yr
        emin  = emean-taumax                ## strain/yr
        azim  = -atan2(x2, ediff) / cov / 2.0e0 ## degrees
        azim  = 90e0+azim
        dexazim = azim+45e0-180e0
        dilat = x1+x3                       ## strain/yr
        sec_inv = sqrt(x1*x1+2e0*x2*x2+x3*x3)
        if params_cov is None:
            staumax, semax, semin, sazim, sdilat, ssec_inv = [None] * 6
        else:
            nv, mv = params_cov.shape
            assert nv == mv and nv == 6
            """
            # cut the part of vcv that holds tau* info
            vcv = params_cov[2:5, 2:5]
            v   = numpy.zeros(shape=(3,1))
            # estimate sigma of tau_max
            v[0,:] = (x1-x3)/4e0/taumax
            v[1,:] = x2/taumax
            v[2,:] = -v[0,:]
            staumax = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
            # estimate sigma of emax
            v[0,:] = 0.5e0*(1+(x1-x3)/2e0/taumax)
            v[1,:] = x2/taumax
            v[2,:] = 0.5e0*(1-(x1-x3)/2e0/taumax)
            semax = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
            # estimate sigma of emin
            v[0,:] = 0.5e0*(1-(x1-x3)/2e0/taumax)
            v[1,:] = -x2/taumax
            v[2,:] = 0.5e0*(1+(x1-x3)/2e0/taumax)
            semin = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
            # estimate sigma of azimuth
            cf = 1e0/((x1-x3)**2e0+4e0*x2**2e0)
            v[0,:] = cf*x2
            v[1,:] = -cf*(x1-x3)
            v[2,:] = -v[0,:]
            sazim = sqrt(numpy.dot(v.T, numpy.dot(vcv, v)))
            # estimate sigma of dilatation
            sdilat = sqrt(vcv[0,0]+vcv[2,2]+2e0*vcv[0,2])
            """
            ## Error propagation for non-linear functions, aka V <- J*VcV*J^T
            ## where J is the Jacobian matrix and VcV the var-covar matrix of
            ## the parameter vector: [Ux, Uy, τx, τxy, τy, ω]
            ## rows of the Jacobian are:
            ## [τ_max, e_max, e_min, Azim, dilatation, sec_inv]
            ## first row is:
            ## [ dτ_max/dUx, dτ_max/dUy, dτ_max/dτx, dτ_max/dτxy, dτ_max/dτy, dτ_max/dω ]
            J = numpy.zeros(shape=(6,6))
            _tmp = ediff/(2e0*taumax)
            J[0, :] = [ 0e0, 0e0, _tmp,           x2/taumax,        -_tmp,          0e0 ]
            J[1, :] = [ 0e0, 0e0, .5e0+_tmp,      x2/taumax,         .5e0-_tmp,     0e0 ]
            J[2, :] = [ 0e0, 0e0, .5e0-_tmp,     -x2/taumax,         .5e0+_tmp,     0e0 ]
            _tmp = ediff*ediff + x2*x2
            J[3, :] = [ 0e0, 0e0, x2/(4e0*_tmp), -ediff/(4e0*_tmp), -x2/(4e0*_tmp), 0e0 ]
            J[4, :] = [ 0e0, 0e0, 1e0,            0e0,               1e0,           0e0 ]
            J[5, :] = [ 0e0, 0e0, x1/sec_inv,     2e0*x2/sec_inv,    x3/sec_inv,    0e0 ]
            Vy = numpy.dot(J, numpy.dot(params_cov, J.T))
            staumax = sqrt(Vy[0,0])
            semax   = sqrt(Vy[1,1])
            semin   = sqrt(Vy[2,2])
            sazim   = sqrt(Vy[3,3])
            sdilat  = sqrt(Vy[4,4])
            ssec_inv= sqrt(Vy[5,5])
        return emean, ediff, \
            taumax, staumax, \
            emax, semax, \
            emin, semin, \
            azim, sazim, \
            dilat, sdilat, \
            sec_inv, ssec_inv

    def info(self):
        return __strain_info__(self.__parameters__)

    def print_details(self, fout, utm_lcm=None):
        """Print Strain Tensor details

            With details, we mean the following parameters:
            'Latitude', 'Longtitude', 'vx+dvx', 'vy+dvy', 'w+dw', 'exx+dexx', \
            'exy+dexy', 'eyy+deyy', 'emax+demax', 'emin+demin', 'shr+dshr', \
            'azi+dazi', 'dilat+ddilat', 'sec. invariant'
            where the units are (in the corresponding order):
            'deg', 'deg', 'mm/yr', 'mm/yr', 'nrad/yr', 'nstrain/yr', \
            'nstrain/yr', 'nstrain/yr', 'nstrain/yr', 'nstrain/yr', \
            'nstrain/yr', 'deg.', 'nstrain/yr', 'nstrain/yr'

            Args:
                fout (output stream): the (already opened) output stream where
                    the details will be written
                utm_zone (int): If given, then the instance's __xcmp__ and 
                    __ycmp__ will be considered UTM Easting and Northing
                    coordinates in the given Zone, and will be transformed to
                    longtitude and latitude before the actual writting takes
                    place.
            Note:
                if the instance's __vcv__ is None (aka we have no var-covar
                matrix), then the sigmas will be printed as '-'
        """
        if utm_lcm:
            cy, cx = [ degrees(c) for c in utm2ell(self.__xcmp__, self.__ycmp__ , None, Ellipsoid("wgs84"), utm_lcm, self.__in_shemisphere__ and utm_lcm>0) ]
        else:
            cx, cy = self.__xcmp__,  self.__ycmp__
        emean, ediff, taumax, staumax, emax, semax, emin, semin, azim, sazim, \
            dilat, sdilat, sec_inv, ssec_inv =  self.cmp_strain(self.__vcv__)
        if self.__vcv__ is not None:
            print('{:9.5f} {:9.5f} {:+7.1f} {:7.1f} {:+7.1f} '\
            '{:7.1f} {:+7.1f} {:7.1f} {:+7.1f} {:7.1f} {:+7.1f} {:7.1f} '\
            '{:+7.1f} {:7.1f} {:+7.1f} {:7.1f} {:+7.1f} {:7.1f} {:+7.1f} '\
            '{:7.1f} {:+7.1f} {:7.1f} {:+7.1f} {:7.1f} {:+7.1f} {:7.1f}'.format(cy, cx, \
            self.value_of('Ux')*1e3, sqrt(self.__vcv__[0,0])*1e3, \
            self.value_of('Uy')*1e3, sqrt(self.__vcv__[1,1])*1e3, \
            self.value_of('omega')*1e9*0.206e0/3.6e0, sqrt(self.__vcv__[5,5])*1e9*0.206e0/3.6e0, \
            self.value_of('taux')*1e9, sqrt(self.__vcv__[2,2])*1e9, \
            self.value_of('tauxy')*1e9, sqrt(self.__vcv__[3,3])*1e9, \
            self.value_of('tauy')*1e9, sqrt(self.__vcv__[4,4])*1e9, \
            emax*1e9, semax*1e9, emin*1e9, semin*1e9, taumax*1e9, \
            staumax*1e9, azim, sazim, dilat*1e9, sdilat*1e9, sec_inv*1e9, ssec_inv*1e9), file=fout)
        else:
            novar = '-'
            print('{:9.5f} {:9.5f} {:+7.1f} {:7s} {:+7.1f} '\
            '{:7s} {:+7.1f} {:7s} {:+7.1f} {:7s} {:+7.1f} {:7s} '\
            '{:+7.1f} {:7s} {:+7.1f} {:7s} {:+7.1f} {:7s} {:+7.1f} '\
            '{:7s} {:+7.1f} {:7s} {:+7.1f} {:7s} {:+7.1f} {:7s}'.format(cy, cx, \
            self.value_of('Ux')*1e3, novar, \
            self.value_of('Uy')*1e3, novar, \
            self.value_of('omega')*1e9*0.206e0/3.6e0, novar, \
            self.value_of('taux')*1e9, novar, \
            self.value_of('tauxy')*1e9, novar, \
            self.value_of('tauy')*1e9, novar, \
            emax*1e9, novar, emin*1e9, novar, taumax*1e9, \
            novar, azim, novar, dilat*1e9, novar, sec_inv*1e9, novar), file=fout)
    
    def print_details_v2(self, fout, utm_lcm=None):
        if utm_lcm:
            #cy, cx = [ degrees(c) for c in utm2ell(self.__xcmp__, self.__ycmp__ , utm_zone) ]
            cy, cx = [ degrees(c) for c in utm2ell(self.__xcmp__, self.__ycmp__ , None, Ellipsoid("wgs84"), utm_lcm, self.__in_shemisphere__ and utm_lcm>0) ]
        else:
            cx, cy = self.__xcmp__,  self.__ycmp__
        emean, ediff, taumax, staumax, emax, semax, emin, semin, azim, sazim, \
            dilat, sdilat, sec_inv, ssec_inv =  self.cmp_strain(self.__vcv__)
        if self.__vcv__ is not None:
            lstr = '%9.5f %9.5f %+7.1f %7.1f %+7.1f '\
            '%7.1f %+7.1f %7.1f %+7.1f %7.1f %+7.1f %7.1f '\
            '%+7.1f %7.1f %+7.1f %7.1f %+7.1f %7.1f %+7.1f '\
            '%7.1f %+7.1f %7.1f %+7.1f %7.1f %+7.1f %7.1f\n' %(cy, cx, \
            self.value_of('Ux')*1e3, sqrt(self.__vcv__[0,0])*1e3, \
            self.value_of('Uy')*1e3, sqrt(self.__vcv__[1,1])*1e3, \
            self.value_of('omega')*1e9*0.206e0/3.6e0, sqrt(self.__vcv__[5,5])*1e9*0.206e0/3.6e0, \
            self.value_of('taux')*1e9, sqrt(self.__vcv__[2,2])*1e9, \
            self.value_of('tauxy')*1e9, sqrt(self.__vcv__[3,3])*1e9, \
            self.value_of('tauy')*1e9, sqrt(self.__vcv__[4,4])*1e9, \
            emax*1e9, semax*1e9, emin*1e9, semin*1e9, taumax*1e9, \
            staumax*1e9, azim, sazim, dilat*1e9, sdilat*1e9, sec_inv*1e9, ssec_inv*1e9)
        else:
            novar = '-'
            lstr = '%9.5f %9.5f %+7.1f %7s %+7.1f '\
            '%7s %+7.1f %7s %+7.1f %7s %+7.1f %7s '\
            '%+7.1f %7s %+7.1f %7s %+7.1f %7s %+7.1f '\
            '%7s %+7.1f %7s %+7.1f %7s %+7.1f %7s\n' %(cy, cx, \
            self.value_of('Ux')*1e3, novar, \
            self.value_of('Uy')*1e3, novar, \
            self.value_of('omega')*1e9*0.206e0/3.6e0, novar, \
            self.value_of('taux')*1e9, novar, \
            self.value_of('tauxy')*1e9, novar, \
            self.value_of('tauy')*1e9, novar, \
            emax*1e9, novar, emin*1e9, novar, taumax*1e9, \
            novar, azim, novar, dilat*1e9, novar, sec_inv*1e9, novar)
        fout.write(lstr)

    def value_of(self, key):
        """Kinda getter.

            This function provides easy access to the instance's attributes,
            specifically the ones in dictionaries. 

            Args:
                Any of the following:
                x returns __xcmp__
                y returns __ycmp__
                any key of the __parameters__ dictionary, or
                any key of the __options__ dictionary
        """
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
        """ Set the x,y values of the instance.

            Args:
                x (float): the value to set __xcmp__ to
                y (float): the value to set __ycmp__ to

        """
        self.__xcmp__ = x
        self.__ycmp__ = y

    def set_to_barycenter(self):
        """ Sets the Strain Tensor (x,y) pair the __stalst__ barycentre
        """
        self.__xcmp__, self.__ycmp__ = barycenter(self.__stalst__)

    def estimate(self):
        """ Estimate fundamental parameters of the Strain Tensor.

            This function will (try to) estimate the Strain Tensor's fundamental
            parameters, aka [Ux, Uy, τx, τxy, τy, ω]. To achieve this, the
            function will perform the following:

            If weighting_function is "shen"
                1. Find the optimal D coefficient
                   Use __options__['d_coef'] if it exists, else find the optimal D
                    coefficient. Assign __options__['d_coef']
                2. Filter stations based on optimal D. The instance's __stalst__
                    will from now on hold only the filtered stations. Assign
                    __stalst__
                3. Compute distance and spatial weights; assign __zweights__ and
                    __lweights__
            4. Compute A and b matrices for LSE (aka Ax = b + ε) and estimate
                x, aka the parameter vector [Ux, Uy, τx, τxy, τy, ω]
            5. Compute a-posteriori std. deviation of the fit and var-covar
                matrix of the estimated parameters. If the latter fails (cause
                it needs an invertion), then the matrix is set to None. Assign
                __vcv__
            6. Assign all fundamental parameters in the __parameters__ dictionary.

            During these steps, the function will assign the following instance
            values:
                self.__sigma0__
                self.__vcv__ (None if stations = 3)
                self.__parameters__['Ux']
                self.__parameters__['Uy']
                self.__parameters__['taux']
                self.__parameters__['tauxy']
                self.__parameters__['tauy']
                self.__parameters__['omega']
            if the weighting_function is shen, then the following members will
            also be assigned:
                self.__zweights__
                self.__lweights__
                self.__stalst__

            Returns:
                numpy.array (6x1): The least squares solution (or if num. of
                    stations is 3 the 'exact' solution) for the Strain.
                    The array contains the estimated values for the parameters:
                    [ Ux, Uy, τx, τxy, τy, ω ]

            Todo:
                When formulating the A and b matrices, we use a sigma0 (aka
                a-priori std. deviation). I think i need this here to compute
                the a-posteriori std. deviation.
        """
        ##  If we are using Shen's weighting sheme, the find D, lweights and
        ##+ zweights.
        if self.__options__['weighting_function'] == 'shen':
            if not self.__options__['d_coef']:
                self.vprint('[DEBUG] Searching for optimal D parameter.')
                if self.__options__['dmin'] >= self.__options__['dmax'] or self.__options__['dstep'] < 0:
                    raise RuntimeError
                lwghts, zwghts, d = self.find_optimal_d()
                self.__stalst__ = self.filter_sta_wrt_distance(d)
            else:
                d = self.__options__['d_coef']
                self.vprint('[DEBUG] Using optimal D parameter {}km.'.format(d))
                self.__stalst__ = self.filter_sta_wrt_distance(d)
                lwghts,_ = self.l_weights()
                zwghts   = self.z_weights()
            self.__zweights__ = zwghts
            self.__lweights__ = lwghts
        ## Formulate the LS matrices A and b (or AW, bW if shen).
        A, b = self.ls_matrices()
        ## Var-Covar matrix
        VcV  = numpy.dot(A.T, A)
        m, n = A.shape
        if m < 6:
            raise RuntimeError('[ERROR] Too few obs to perform LS.')
        elif m == 6:
            self.vprint('[DEBUG] Only 3 stations available; computing NOT estimating strain.')
            self.__vcv__ = None
        ##  Note: To silence warning in versions > 1.14.0, use a third argument,
        ##+ rcond=None; see https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
        estim, res, rank, sing_vals = numpy.linalg.lstsq(A, b)
        # Parameter variance-covariance matrix
        if m > 6:
            try:
                ##  A-posteriori std. deviation. res[0] is the sum of residuals;
                ##+ squared Euclidean 2-norm for each column in b - a*x, aka
                ##+ u^T * P * u,
                ##+ or (b - A*estim)^T * (b - A*estim)
                sigma0_post = float(res[0])
                self.__sigma0__ = sqrt(sigma0_post / (float(m) - 6e0))
                bvar = linalg.inv(VcV) * (sigma0_post/float(m-n))
                ## A-posteriri VcV matrix = σ0^2 * (A^T P A)^-1
                self.__vcv__ = bvar
            except:
                self.vprint('[DEBUG] Cannot compute var-covar matrix! Probably singular.')
                self.__vcv__ = None
        self.__parameters__['Ux']    = float(estim[0])
        self.__parameters__['Uy']    = float(estim[1])
        self.__parameters__['taux']  = float(estim[2])
        self.__parameters__['tauxy'] = float(estim[3])
        self.__parameters__['tauy']  = float(estim[4])
        self.__parameters__['omega'] = float(estim[5])
        return estim
