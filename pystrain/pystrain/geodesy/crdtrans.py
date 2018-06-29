#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
from pystrain.geodesy.ellipsoid import Ellipsoid

def top2daz(north, east, up):
    """Compute azimouth, zenith and distance from a topocentric vector.

        Given a topocentric vector (aka north, east and up components), compute
        the azimouth, zenith angle and distance between the two points.

        Args:
            north (float): the north component (in meters)
            east (float) : the east component (in meters)
            up (float)   : the up component (in meters)

        Returns:
            tuple (floats): a tuple of three floats is returned, as:
                            [distance, azimouth, zenith], where distance is
                            in meters, and azimouth and zenith are in radians.
    """
    distance = math.sqrt(north*north + east*east + up*up)
    a        = math.atan2(east, north) % (math.pi*2e0) # normalized [0, 2pi]
    zenith   = math.acos(up/distance);
    return distance, a, zenith

def car2top(xi, yi, zi, xj, yj, zj, ell=Ellipsoid("wgs84")):
    """Cartesian to topocentric vector.

        Transform a vector expressed in cartesian coordinates to the topocentric,
        local system around point i (i.e. North(i), East(i), Up(i)).
        To perform the conversion we will need an ellipsoid. The default is
        'wgs84'.

        Args:
            xi (float): x-component of reference point (m).
            yi (float): y-component of reference point (m).
            zi (float): z-component of reference point (m).
            xj (float): x-component of end point (m).
            yj (float): y-component of end point (m).
            zj (float): z-component of end point (m).

        Returns:
            tuple (float): a 3-element float tuple, as [north, east, up] in
                           meters.

        Note:
            The vector transformed is [xj-xi, yj-yi, zj-zi]
    """
    # Cartesian to ellipsoidal for reference point.
    phi_i, lamda_i, h_i = car2ell(xi, yi, zi, ell)

    # Trigonometric numbers.
    cosf = math.cos(phi_i)
    cosl = math.cos(lamda_i)
    sinf = math.sin(phi_i)
    sinl = math.sin(lamda_i)

    # Cartesian vector.
    dx = xj - xi
    dy = yj - yi
    dz = zj - zi

    # Topocentric vector.
    north = - sinf * cosl * dx - sinf * sinl * dy + cosf * dz
    east  = - sinl * dx        + cosl * dy
    up = cosf * cosl * dx + cosf * sinl * dy + sinf * dz

    return north, east, up

def ell2car(phi, lamda, h, ell=Ellipsoid("wgs84")):
    """Ellipsoidal to cartesian coordinates.

        Convert Ellipsoidal coordinates (aka longtitude, latitude, height) to
        cartesian. Default ellipsoid is wgs84.

        Args:
            phi (float)    : the latitude (in radians).
            lamda (float)  : the longtitude (in radians).
            h (float)      : the height (in meters)
            ell (Ellipsoid): the ellipsoid of choice.
        
        Returns: 
            tuple (float): a 3-float tuple, as [x, y, z] in meters
    """
    # Eccentricity squared.
    e2 = ell.eccentricity_squared()

    # Trigonometric numbers.
    sinf = math.sin(phi)
    cosf = math.cos(phi)
    sinl = math.sin(lamda)
    cosl = math.cos(lamda)

    # Radius of curvature in the prime vertical.
    N = ell.N(phi)

    # Compute geocentric rectangular coordinates.
    x = (N+h) * cosf * cosl;
    y = (N+h) * cosf * sinl;
    z = ((1.0e0-e2) * N + h) * sinf;

    # Finished.
    return x, y, z

def car2ell(x, y, z, ell=Ellipsoid("wgs84")):
    """Cartesian to ellipsoidal coordinates.
    
        Cartesian to ellipsoidal coordinates (aka latitude, longtitude and
        height). Reference: Fukushima, T., "Transformation from Cartesian to 
        geodetic coordinates accelerated by Halley's method",
        J. Geodesy (2006), 79(12): 689-693

        Args:
            x (float): cartesian x component (meters).
            y (float): cartesian y component (meters).
            z (float): cartesian z component (meters).
            ell      : ellipsoid of choice.

        Returns:
            tuple (float): the ellipsoidal coordinates in a 3-element tuple, as
                           [latitude, longtitude, height]. Longtitude and
                           latitude are returned in radians, while height is in
                           meters.
    """
    a = ell.a
    f = ell.f
    # Functions of ellipsoid parameters.
    aeps2 = a*a*1e-32
    e2    = (2.0e0-f)*f
    e4t   = e2*e2*1.5e0
    ep2   = 1.0e0-e2
    ep    = math.sqrt(ep2)
    aep   = a*ep

    # Compute distance from polar axis squared.
    p2 = x*x + y*y

    # Compute longitude lamda.
    if  (p2 != 0):
        lamda = math.atan2(y,x);
    else:
        lamda = .0e0;

    # Ensure that Z-coordinate is unsigned.
    absz = abs(z)

    if (p2 > aeps2): # Continue unless at the poles
        # Compute distance from polar axis.
        p = math.sqrt(p2)
        # Normalize.
        s0  = absz/a
        pn  = p/a
        zp  = ep*s0
        # Prepare Newton correction factors.
        c0  = ep*pn
        c02 = c0*c0
        c03 = c02*c0
        s02 = s0*s0
        s03 = s02*s0
        a02 = c02+s02
        a0  = math.sqrt(a02)
        a03 = a02*a0
        d0  = zp*a03 + e2*s03
        f0  = pn*a03 - e2*c03
        # Prepare Halley correction factor.
        b0  = e4t*s02*c02*pn*(a0-ep)
        s1  = d0*f0 - b0*s0
        cp  = ep*(f0*f0-b0*c0)
        # Evaluate latitude and height.
        phi = math.atan(s1/cp);
        s12 = s1*s1
        cp2 = cp*cp
        h = (p*cp+absz*s1-a*math.sqrt(ep2*s12+cp2))/math.sqrt(s12+cp2)

    else: # Special case: pole.
        phi = math.pi / 2e0
        h   = absz - aep

    # Restore sign of latitude.
    if (z < 0.e0):
        phi = -phi

    # Finished.
    return phi, lamda, h

if __name__ == "__main__":
    dyng_xyz = [4595220.002e0, 2039434.077e0, 3912625.997e0]
    lat, lon, hgt = car2ell(dyng_xyz[0], dyng_xyz[1], dyng_xyz[2])
    dx, dy, dz = ell2car(lat, lon, hgt)
    print('From cartesian to ellipsoidal and back, diffs are:')
    print('Δx = {}\nΔy = {}\nΔz = {}'.format(abs(dyng_xyz[0]-dx),
        abs(dyng_xyz[1]-dy), abs(dyng_xyz[2]-dz)))
    assert abs(dyng_xyz[0]-dx) < 1e-7 and \
        abs(dyng_xyz[1]-dy) < 1e-7 and \
        abs(dyng_xyz[2]-dz) < 1e-7

    dyng2xyz = [ c+5e0 for c in dyng_xyz ]
    lat2, lon2, hgt2 = car2ell(dyng2xyz[0], dyng2xyz[1], dyng2xyz[2])
    n1, e1, u1 = car2top(dyng_xyz[0], dyng_xyz[1], dyng_xyz[2],
        dyng2xyz[0], dyng2xyz[1], dyng2xyz[2])
    print('Cartesian to Topocentric')
    print('Δn = {}\nΔe = {}\nΔu = {}'.format(n1,e1,u1))
    dr = math.sqrt(3e0 * math.pow(5e0,2))
    assert abs(dr - math.sqrt(n1*n1 + e1*e1 +u1*u1)) < 1e-5
    distance, a, zenith = top2daz(n1, e1, u1)
    assert abs(dr-distance) < 1e-5
