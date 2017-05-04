#! /usr/bin/python

from math import floor, degrees, radians, pi, sin, cos, tan
from pystrain.geodesy.ellipsoid import Ellipsoid

def ell2utm(lat, lon, ell=Ellipsoid("wgs84"), lcm=None):
    f  = ell.f
    a  = ell.a
    e2 = ell.eccentricity_squared()

    Zone = floor(degrees(lon)/6)+31
    Zone = Zone + int(Zone<=0)*60 - int(Zone>60)*60
    lcm  = radians(Zone*6-183)

    ko = 0.9996   # Scale factor
    if lat > 0:
        No = 0    # False northing (north)
    else:
        No = 1e7  # False northing (south)
    Eo = 500000   # False easting

    lam = lon-lcm;
    if lam >= pi: lam = lam - pi*2

    RN = ell.N(lat)
    RM = ell.M(lat)

    coslat = cos(lat)
    sinlat = sin(lat)
    h2     = e2*coslat*coslat/(1-e2)
    t      = tan(lat)
    n      = f/(2-f)

    # powers of various values
    n2   = pow(n,2)
    n3   = n2*n
    n4   = n3*n
    t2   = pow(t,2)
    t3   = t2*t
    t4   = t3*t
    t6   = t4*t2
    h22  = pow(h2,2)
    h23  = h22*h2
    h24  = h23*h2

    A0 = 1+n2/4.0+n4/64.0
    A2 = 3.0/2.0*(n-n3/8)
    A4 = 15.0/16.0*(n2-n4/4)
    A6 = 35.0/48.0*n3
    A8 = 315.0/512.0*n4
    S  = a/(1+n)*(A0*lat-A2*sin(2*lat)+A4*sin(4*lat)-A6*sin(6*lat)+A8*sin(8*lat))

    E1   = lam*coslat
    E2   = pow(lam,3)*pow(coslat,3)/6*(1-t2+h2)
    E3   = pow(lam,5)*pow(coslat,5)/120*(5-18*t2+t4+14*h2-58*t2*h2+
               13*h22+4*h23-64*t2*h22-24*t2*h23)
    E4   = pow(lam,7)*pow(coslat,7)/5040*(61-479*t2+179*t4-t4*t2)
    E    = Eo+ko*RN*(E1+E2+E3+E4)

    N1 = S/RN
    N2 = pow(lam,2)/2*sinlat*coslat;
    N3 = pow(lam,4)/24*sinlat*pow(coslat,3)*(5-t2+9*h2+4*h22)
    N4 = pow(lam,6)/720*sinlat*pow(coslat,5)*(61-58*t2+t4+
            270*h2-330*t2*h2+445*h22+324*h23-680*t2*h22+
            88*h24-600*t2*h23-192*t2*h24)
    N5 = pow(lam,8)/40320*sinlat*pow(coslat,7)*(1385-311*t2+543*t4-t6)
    N = No+ko*RN*(N1+N2+N3+N4+N5)

    return N, E, Zone, lcm

if __name__ == "__main__":
    # lats = [39.010444, -12.0464, -37.8136, 38.7223]
    lats = [0.6808606982218858, -0.2102493430122449, -0.6599718220321278, 0.6758316289450003]
    # lons = [22.969694, -77.0428, 144.9631, -9.1393]
    lons = [0.40089679623260527, -1.3446505249554874, 2.530083388897792, -0.15951087632751776]
    # ellipsoid grs80
    ell = Ellipsoid('grs80')
    # octave results
    octave = [
        [4319781.22749958, 670541.267192695, 34, 0.366519142918809],
        [8667487.89699582, 277617.453174038, 18, -1.30899693899575],
        [5812911.69963288, 320704.446318808, 55, 2.56563400043166],
        [4285969.85875752, 487890.758080714, 29, -0.157079632679490]]
    for i in range(0, len(lats)):
        n, e, z, l = ell2utm(lats[i], lons[i], ell)
        print '{} {} {} {}'.format(n, e, z, l)
        print 'Octave diffs={} {} {} {}'.format(abs(n-octave[i][0]), abs(e-octave[i][1]), abs(z-octave[i][2]), abs(l-octave[i][3]))
