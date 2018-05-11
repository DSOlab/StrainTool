""" This module defines a Station class. A station is supposed to represent a
    point on the globe. It has coordinates (usually defined as longtitude and
    latitude), a name, (tectonic) velocities and respective standard deviations
    (in east and north components), a correlation coefficient between East and
    North velocity components and a time-span.
    This class is designed to assist the estimation of strain tensors; hence,
    only attributes that could help with this are considered.
"""

from math import sqrt, radians

# Any Station instance, can have any (or all) of these attributes
station_member_names = ['name', 'lat', 'lon', 've', 'vn', 'se', 'sn', 'rho', 't']

class Station:

    def __init__(self, *args, **kargs):
        """ Station constructor; construction can be performed:
            1 from an input string of type:
                "name lon lat Ve Vn Se Sn RHO T"
                where lon and lat are in decimal degrees and velocity components
                are in mm/year.
            2 given any of the (above mentioned) instance members/attributes.

            e.g. s = Station("akyr +24.91260690 +34.98083160 0.00871244 -0.0151236 0.00136367 0.000278371 0.5  2.5")
                 s = Station(name="akyr")
                 s = Station(name="akyr", lat=34.98083160, ve=-0.0151236)
        """
        self.set_none()

        if len(args) is not 0:
            self.init_from_ascii_line(args[0])

        if len(kargs) is not 0:
            for key, val in kargs.items():
                if key in station_member_names:
                    setattr(self, key, val)

    def init_from_ascii_line(self, input_line):
        l = input_line.split()
        try:
            self.name = l[0]
            self.lon  = radians(float(l[1]))
            self.lat  = radians(float(l[2]))
            self.ve   = float(l[3]) / 1000e0
            self.vn   = float(l[4]) / 1000e0
            self.se   = float(l[5]) / 1000e0
            self.sn   = float(l[6]) / 1000e0
            self.rho  = float(l[7])
            self.t    = float(l[8])
        except:
            print '[DEBUG] Invalid Station instance constrution.'
            print '[DEBUG] Input line \"{}\"'.format(input_line.strip())
            raise RuntimeError

    def set_none(self):
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
        dlon = sta.lon - self.lon
        dlat = sta.lat - self.lat
        return dlon, dlat, sqrt(dlat*dlat + dlon*dlon)
