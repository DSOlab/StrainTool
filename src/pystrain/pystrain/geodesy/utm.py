#! /usr/bin/python

from math import floor, degrees, radians, pi, sin, cos, tan
from pystrain.geodesy.ellipsoid import Ellipsoid

def utm2ell(E, N, zone, ell=Ellipsoid("wgs84"), lcm=None):
    if not lcm:
        lcm = radians(abs(zone)*6-183)

    f   = ell.f
    a   = ell.a
    e2  = ell.eccentricity_squared()
    e22 = e2*e2
    e23 = e22*e2
    e24 = e23*e2

    No = 0           # False northing (north)
    if zone < 0:
        No = 1e7     # False northing (south)
    Eo = 500000      # False easting
    N = N-No
    E = E-Eo
    Zone = abs(zone) # Remove negative zone indicator for southern hemisphere
    ko = 0.9996      # UTM scale factor
    lat1 = N/ko/a
    dlat = 1
    while abs(dlat) > 1e-12:
        A0=1-(e2/4)-(e22*3/64.0)-(e23*5/256.0)-(e24*175/16384.0)
        A2=(3/8.0)*( e2+(e22/4.0)+(e23*15/128.0)-(e24*455/4096.0) )
        A4=(15/256.0)*( e22+(e23*3/4.0)-(e24*77/128.0) )
        A6=(35/3072.0)*( e23-(e24*41/32.0) )
        A8=-(315/131072.0)*e24
        f1=a*( A0*lat1-A2*sin(2*lat1)+A4*sin(4*lat1)-A6*sin(6*lat1)+A8*sin(8*lat1) )-N/ko
        f2=a*( A0-2*A2*cos(2*lat1)+4*A4*cos(4*lat1)-6*A6*cos(6*lat1)+8*A8*cos(8*lat1) )
        dlat=-f1/f2
        lat1=lat1+dlat
    RN = ell.N(lat)
    RM = ell.M(lat)
    h2 = e2*pow(cos(lat1),2)/(1-e2)
    t  = tan(lat1)
    t2 = pow(t,2)
    t4 = t2*t2
    t6 = t4*t2

    E0 = E/ko/RN
    E1 = E0
    E2 = pow(E0,3)/6.*(1+2*t2+h2)
    E3 = pow(E0,5)/120.*(5+6*h2+28*t2-3*h2^2+8*t2*h2+24*t4-4*h2^3+4*t2*h2^2+24*t2*h2^3)
    E4 = pow(E0,7)/5040.*(61 + 662*t2 + 1320*t4 + 720*t6)
    lon = sec(lat1)*(E1-E2+E3-E4)+lcm

    E0=E/ko
    N1=(t.*E0.^2)./(2*RM.*RN)
    N2=(t.*E0.^4)./(24*RM.*RN.^3).*(5+3.*t.^2+h2-4*h2.^2-9*h2.*t.^2)
    N3=(t.*E0.^6)./(720*RM.*RN.^5).*(61-90*t.^2+46*h2+45*t.^4-252*t.^2.*h2- ...
            5*h2.^2+100*h2.^3-66*t.^2.*h2.^2-90*t.^4.*h2+88*h2.^4+225*t.^4.*h2.^2+ ...
            84*t.^2.*h2.^3 - 192*t.^2.*h2.^4)
    N4=(t.*E0.^8)./(40320*RM.*RN.^7).*(1385+3633*t.^2+4095*t.^4+1575*t.^6)
    lat=lat1-N1+N2-N3+N4

def ell2utm(lat, lon, ell=Ellipsoid("wgs84"), zone=None, lcm=None):
    """
        All input arguments in radians.
        TODO zone should be given in degrees
    """
    f  = ell.f
    a  = ell.a
    e2 = ell.eccentricity_squared()

    if zone:
        Zone = zone
    else:
        Zone = floor(degrees(lon)/6)+31
        Zone = Zone + int(Zone<=0)*60 - int(Zone>60)*60
    lcm  = radians(Zone*6-183)

    ko = 0.9996   # Scale factor
    if lat > 0:
        No = 0    # False northing (north)
    else:
        No = 1e7  # False northing (south)
    Eo = 500000   # False easting

    lam = lon-lcm
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
