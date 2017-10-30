#! /bin/python2.7

from math import cos, sin, sqrt

##  A dictionary holding standard reference ellipsoids. The keys are the
##+ ellipsoid names, and the respective values are the defining geometric
##+ parameters, i.e. a (semi-major axis) and f (flattening).
ref_ell_dict = {
    "grs80": [6378137.0e0, 1.0e00/298.257222101e0],
    "wgs84": [6378137.0e0, 1.0e00/298.257223563e0],
    "pz90" : [6378135.0e0, 1.0e00/298.257839303e0]
}

class Ellipsoid:
    """
        A class to represent reference ellipsoids.
    """
    def __init__(self, name, a=None, f=None):
        """
            Constructor for Ellipsoid. You can either pass in a standard name
            (i.e. included in the ref_ell_dict dictionary), or construct a
            custom one, using the 'a' and 'f' parameters.
            If you do use the 'a' and 'f' parameters, make sure that the 'name'
            is not a standard one.
        """
        # a and f can be both None or both not None
        if bool(a) + bool(f) == 1:
            raise RuntimeError
        self.name = name
        # if a and f are None, the user (probably) wants a standard ellipsoid
        if not a:
            if not ref_ell_dict.has_key(name.lower()):
                raise LookupError
            self.a, self.f = ref_ell_dict[name.lower()]
        # user has supplied a and f vals
        else:
            if ref_ell_dict.has_key(name.lower()):
                raise RuntimeError
            self.a = a
            self.f = f

    def eccentricity_squared(self):
        """ Computed and return the ellipsoid's squared eccentricity.
        """
        return ( 2.0e0 - self.f ) * self.f;

    def semi_minor(self):
        """ Computed and return the ellipsoid's semi-minor axis.
        """
        return self.a * (1.0e0 - self.f)

    def N(self, lat):
        """  Compute the normal radius of curvature at a given latitude.
            'lat' parameter should be in radians.
        """
        cosf  = cos(lat)
        sinf  = sin(lat)
        acosf = self.a * cosf
        bsinf = sinf * self.semi_minor()
        den   = sqrt(acosf*acosf + bsinf*bsinf)
        return (self.a * self.a) / den;

    def M(self, lat):
        """  Compute the meridional radii of curvature at a given latitude.
            'lat' parameter should be in radians.
        """
        a     = self.a
        b     = self.semi_minor()
        cosf  = cos(lat)
        sinf  = sin(lat)
        acosf = a * cosf
        bsinf = b * sinf
        tmpd  = acosf*acosf + bsinf*bsinf
        return ( (a*b)/tmpd ) * ( (a*b)/sqrt(tmpd) )

    def __getattr__(self, name):
        """ For ease of use, the instances has the following attributes:
                * e2   -> eccentricity squared
                * b    -> semi-minor axis
                * finv -> inverse flattening
            I.e. the user can legaly write:
            ell1 = Ellipsoid("grs80")
            ell1.b
            Of course, these attributes cannot be assigned to.
        """
        # print 'Attribute not found! Got in getattr().'
        if name == "e2"  : return self.eccentricity_squared()
        if name == "b"   : return self.semi_minor()
        if name == "finv": return 1.0e0/self.f
        raise AttributeError

if __name__ == "__main__":
    ell1 = Ellipsoid("grs80")
    ell2 = Ellipsoid("GRS80")
    ell3 = Ellipsoid("mine", 1.2, 2.3)
    try:
        # this should throw
        ell4 = Ellipsoid("mine", 1.2)
    except:
        print 'Caught exception for ell4.'
    try:
        ell5 = Ellipsoid("pz90", 1.2, 4.5)
    except:
        print 'Caught exception for ell5.'
    print 'Normal radius of curvature at ~ Athens {:15.3f}'.format(ell1.N(0.6628663536338242))
    print 'Meridional radius of curvature at ~ Athens {:15.3f}'.format(ell1.M(0.6628663536338242))
    print 'grs80 geometric parameters:'
    print '\tSemi-major           {:10.3f}'.format(ell1.a)
    print '\tSemi-minor           {:10.3f}'.format(ell1.b)
    print '\tFLattening           {:10.3f}'.format(ell1.f)
    print '\tInv. Flattening      {:10.3f}'.format(ell1.finv)
    print '\tEccentricity Squared {:10.3f}'.format(ell1.e2)
