# -*- coding: utf-8 -*-

from math import sqrt, radians, sin, cos, atan2, pi

# Any Station instance, can have any (or all) of these attributes
station_member_names = ['name', 'lat', 'lon', 've', 'vn', 'se', 'sn', 'rho', 't']

class Station:
    '''A simple Station class.

        This module defines a Station class. A station is supposed to represent a
        point on the globe. It has coordinates (usually defined as longtitude and
        latitude), a name, (tectonic) velocities and respective standard deviations
        (in east and north components), a correlation coefficient between East and
        North velocity components and a time-span.
        This class is designed to assist the estimation of strain tensors; hence,
        only attributes that could help with this are considered.

        Attributes:
            name (str) : the name of the station
            lon (float): longtitude of the station (radians). In case the station
                         coordinates are transformed to easting and northing
                         (aka to projection coordinates), this component will
                         hold the Easting.
            lat (float): latitude of the station (radians). In case the station
                         coordinates are transformed to easting and northing
                         (aka to projection coordinates), this component will
                         hold the Northing.

            ve (float) : velocity of the east component, in meters/year
            vn (float) : velocity of the north component, in meters/year
            se (float) : std. deviation of the east velocity component, in
                         meters/year
            sn (float) : std. deviation of the north velocity component, in
                         meters/year
            rho (float): correlation coefficient between East and North velocity
                         components
            t (float)  : time-span in decimal years
            
    '''

    def __init__(self, *args, **kargs):
        '''Station constructor.
        
            Station constructor; construction can be performed:
                #. from an input string of type:
                    "name lon lat Ve Vn Se Sn RHO T"
                    where lon and lat are in decimal degrees and velocity components
                    are in mm/year.
                #. given any of the (above mentioned) instance members/attributes.

            e.g. s = Station("akyr +24.91260690 +34.98083160 8.71244 -15.1236 0.00136367 0.000278371 0.5  2.5")
                 s = Station(name="akyr")
                 s = Station(name="akyr", lat=34.98083160, ve=-0.0151236)

            Args:
                *args (str): if provided, then it is supposed to be a station
                             string ("name lon lat Ve Vn Se Sn RHO T") and the
                             function will try to resolve it and assign member
                             values.
                **kargs:     any named member variable, aka one of:
                    * name
                    * lon
                    * lat
                    * ve
                    * vn
                    * se
                    * sn
                    * rho
                    * t
        '''
        self.set_none()

        if len(args) is not 0:
            self.init_from_ascii_line(args[0])

        if len(kargs) is not 0:
            for key, val in kargs.items():
                if key in station_member_names:
                    setattr(self, key, val)

    def init_from_ascii_line(self, input_line):
        '''Assignment from string.

            This function will initialize all member values of a station
            instance, given a (string) line of type:
            "name lon lat Ve Vn Se Sn RHO T"
            where lon and lat are in decimal degrees and velocity components
            and sigmas are in mm/year.

            Args:
                input_line (str): a string (line) of type:
                                  "name lon lat Ve Vn Se Sn RHO T"

            Raises:
                RuntimeError: if the input line (string) cannot be resolved
        '''
        l = input_line.split()
        try:
            self.name = l[0]
            self.lon  = radians(float(l[1]))
            self.lat  = radians(float(l[2]))
            self.ve   = float(l[3]) / 1e3
            self.vn   = float(l[4]) / 1e3
            self.se   = float(l[5]) / 1e3
            self.sn   = float(l[6]) / 1e3
            self.rho  = float(l[7]) / 1e3
            self.t    = float(l[8])
        except:
            print('[DEBUG] Invalid Station instance constrution.')
            print('[DEBUG] Input line \"{}\"'.format(input_line.strip()))
            raise RuntimeError

    def set_none(self):
        '''Set to None.

            Set all instance member values to None.
        '''
        self.name = None                                                        
        self.lon  = None
        self.lat  = None
        self.ve   = None
        self.vn   = None
        self.se   = None
        self.sn   = None
        self.rho  = None
        self.t    = None

    def distance_from(self, sta):
        '''Distance to another station.

            Compute the distance of the instance to another instance (of type
            Station). The algorithm is just an Eucledian norm, so the actual
            station components must already have been transformed to cartesian.

            Args:
                sta (Station): a station instance

            Returns:
                tuple (float, float, float): a 3-float tuple, where the elements
                are
                    #. dlon 
                    #. dlat
                    #. dr
                If the calling station has index i and the station passed in
                has index j, then the returned values are computed as
                    * δlon = lon_j - lon_i
                    * δlat = lat_j - lat_i
                    * δr   = sqrt{δlon**2 + δlat**2}
              
            Warning:
                For the function to return valid results, the station coordinate
                component must not be in ellipsoidal coordinates; the station
                coordinates should have already been transformed to cartesian
                before calling this function. The function will treat the "lon"
                attribute as "x" or "Easting" and the "lat" component as "y" or
                "Northing".

        '''
        dlon = sta.lon - self.lon
        dlat = sta.lat - self.lat
        return dlon, dlat, sqrt(dlat*dlat + dlon*dlon)

    def haversine_distance(self, sta, R=6371008.8e0):
        """Distance between two point on sphere.

            This function will calculate the distance on a spherical earth
            (ignoring ellipsoidal effects), via the ‘haversine’ formula 
            (great-circle distance between two points). Note that  using a 
            spherical model gives errors typically up to 0.3%. 
            If the calling station has index i and the station passed in
            has index j, then the returned value is computed as
            * δr = sta.lon - self.lon

            Args:
                sta (Station): a station instance
                R (float): mean Earth's radius in meters

            Returns:
                tuple (float, float): First element is distance in m and second
                    is the error, computed as 0.3%.

            Warning:
                It is assumed that the stations's (both calling instance and sta)
                .lon and .lat components are in radians.

        """
        if self.lon < 0e0 or sta.lon < 0e0:
            selflon, stalon = self.lon + 2e0*pi, sta.lon  + 2e0*pi
        else:
            selflon, stalon = self.lon, sta.lon
        Df = sta.lat - self.lat
        Dl = stalon - selflon
        sinDf2 = sin(Df/2e0)
        cosDf2 = cos(Df/2e0)
        a  = sinDf2*sinDf2 + cos(sta.lat)*cos(self.lat)*sinDf2*sinDf2
        c  = 2e0*atan2(sqrt(a), sqrt(1e0-a))
        dr = R*c
        ddr = (dr * 0.3e0) / 100e0
        return dr, ddr
