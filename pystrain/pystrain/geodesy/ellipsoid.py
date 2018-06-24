#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from math import cos, sin, sqrt

##  A dictionary holding standard reference ellipsoids. The keys are the
##+ ellipsoid names, and the respective values are the defining geometric
##+ parameters, i.e. a (semi-major axis) and f (flattening).
ref_ell_dict = {
    "grs80": [6378137e0, 1e0/298.257222101e0],
    "wgs84": [6378137e0, 1e0/298.257223563e0],
    "pz90" : [6378135e0, 1e0/298.257839303e0]
}

class Ellipsoid:
    """A class to represent reference ellipsoids.

        This class constructs reference ellipsoids, see
        https://en.wikipedia.org/wiki/Reference_ellipsoid

        Attributes:
            name (str): the name of the ellipsoid.
            a (float) : the ellipsoid's semi-major axis (meters).
            f (float) : the ellipsoid's flattening.
    """

    def __init__(self, name, a=None, f=None):
        """Ellipsoid constructor.

            Constructor for Ellipsoid. You can either pass in a standard name
            (i.e. included in the ref_ell_dict dictionary), or construct a
            custom one, using the 'a' and 'f' parameters.
            If you do use the 'a' and 'f' parameters, make sure that the 'name'
            is not a standard one.

            Args:
                name (str): the ellipsoid's name; if a and f are not passed in,
                            then the function will try to match this parameter
                            to a reference ellipsoid in the ref_ell_dict
                            dictionary. If a and f are passed in, then a custom
                            ellipsoid is constructed.
                a (float) : The semimajor axis in meters; only use this parameter
                            to construct a custom ellipsoid. If this parameter
                            is used, then the user must also specify the f parameter.
                b (float) :The flattening of the ellipsoid; only use this parameter
                            to construct a custom ellipsoid. If this parameter
                            is used, then the user must also specify the a parameter.

            Raises:
                LookupError: if the parameters a and f are not both specified, or
                             both not specified.
                RuntimeError: if the user has specified a name for the ellipsoid
                              that is not in the ref_ell_dict dictionary (without
                              providing the a and f parameters).
        """
        # a and f can be both None or both not None
        if bool(a) + bool(f) == 1:
            raise RuntimeError
        self.name = name
        # if a and f are None, the user (probably) wants a standard ellipsoid
        if not a:
            #if not ref_ell_dict.has_key(name.lower()): !! not Python3 compatible
            if name.lower() not in ref_ell_dict:
                raise LookupError
            self.a, self.f = ref_ell_dict[name.lower()]
        # user has supplied a and f vals
        else:
            #if ref_ell_dict.has_key(name.lower()):
            if name.lower() in ref_ell_dict:
                raise RuntimeError
            self.a = a
            self.f = f

    def eccentricity_squared(self):
        """Compute and return the ellipsoid's squared eccentricity.

            Returns:
                float: the squared eccentricity.
        """
        return (2e0-self.f)*self.f;

    def semi_minor(self):
        """Compute and return the ellipsoid's semi-minor axis.

            Returns:
                float: the ellipsoid's semi-minor axis (meters).
        """
        return self.a*(1e0-self.f)

    def N(self, lat):
        """Normal radius of curvature.
        
            Given a latitude on the ellipsoid, compute the normal radius of
            curvature (on the parallel).

            Args:
                lat (float): the latitude in radians.

            Returns:
                float: normal radius of curvature at given latitude (meters).
        """
        cosf  = cos(lat)
        sinf  = sin(lat)
        acosf = self.a * cosf
        bsinf = sinf * self.semi_minor()
        den   = sqrt(acosf*acosf + bsinf*bsinf)
        return (self.a * self.a) / den;

    def M(self, lat):
        """Meridional radii of curvature.
        
            Compute the meridional radii of curvature at a given latitude
            (on the parallel).

            Args:
                lat (float): the latitude in radians.

            Returns:
                float: the meridional radii of curvature at the given latitude
                       (meters).
        """
        a     = self.a
        b     = self.semi_minor()
        cosf  = cos(lat)
        sinf  = sin(lat)
        acosf = a * cosf
        bsinf = b * sinf
        tmpd  = acosf*acosf + bsinf*bsinf
        return ((a*b)/tmpd) * ((a*b)/sqrt(tmpd))

    def __getattr__(self, name):
        """(Attribute) getter.
        
            For ease of use, the instances has the following attributes:
                * e2   -> eccentricity squared
                * b    -> semi-minor axis
                * finv -> inverse flattening
            I.e. the user can legaly write:
            ell1 = Ellipsoid("grs80")
            ell1.b
            Of course, these attributes cannot be assigned to.

            Args:
                name (str): any (valid) attribute we want to get.

            Returns:
                
        """
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
        print('Caught exception for ell4.')
    try:
        ell5 = Ellipsoid("pz90", 1.2, 4.5)
    except:
        print('Caught exception for ell5.')
    print('Normal radius of curvature at ~ Athens {:15.3f}'.format(ell1.N(0.6628663536338242)))
    print('Meridional radius of curvature at ~ Athens {:15.3f}'.format(ell1.M(0.6628663536338242)))
    print('grs80 geometric parameters:')
    print('\tSemi-major           {:10.3f}'.format(ell1.a))
    print('\tSemi-minor           {:10.3f}'.format(ell1.b))
    print('\tFLattening           {:10.3f}'.format(ell1.f))
    print('\tInv. Flattening      {:10.3f}'.format(ell1.finv))
    print('\tEccentricity Squared {:10.3f}'.format(ell1.e2))
