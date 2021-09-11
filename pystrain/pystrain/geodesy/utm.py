#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from pystrain.geodesy.ellipsoid import Ellipsoid
import math

MAX_UTM_ITERATIONS = 100

def dd2dms(dd):
    ''' Decimal degrees to hexicondal degrees.
        
        Args:
            dd (float): angle in decimal degrees

        Returns:
            tuple (int, int, float): The elements are:
                # 0 -> integer degrees
                # 1 -> integer minutes
                # 2 -> float seconds
    '''
    dd1  = abs(float(dd))
    cdeg = int(dd1)
    minsec = dd1 - cdeg
    cmin = int(minsec * 60)
    csec = (minsec % 60) / float(3600)
    if dd < 0e0: cdeg = cdeg * -1
    return cdeg,cmin,csec

def zone2lcm(zone_number):
  return (zone_number-1)*6-180+3

def utm2ell(E, N, zone=None, ell=Ellipsoid("wgs84"), lcm=None, southern_hemisphere=False):
    '''UTM to ellipsoidal coordinates.

        Convert UTM coordinates (i.e. Easting and Northing) to ellipsoidal
        coordinates (actualy latitude and longtitude).

        Args:
            E (float): Easting in meters
            N (float): Northing in meters
            zone (int): the zone in degrees (ignored if lcm set)
            ell (Ellipsoid): the ellipsoid of choice
            lcm: non-standard central meridian (radians)

        Returns:
            tuple (float, float): first is latitude and second is longtitude,
                                  both in degrees.
    '''
    if zone == lcm == None:
      raise RuntimeError('[ERROR] utm2ell:: need to specify at least zone or non-standard central meridian')

    if not lcm:
        lcm = math.radians(abs(zone)*6-183)

    if southern_hemisphere:
      N -= 10000000e0

    f   = ell.f
    a   = ell.a
    e2  = ell.eccentricity_squared()
    e22 = e2*e2
    e23 = e22*e2
    e24 = e23*e2

    No = 1e7 if lcm < 0 else 0e0 # False northing (north)
    Eo   = 500000    # False easting
    N    = N-No
    E    = E-Eo
    # Zone = abs(zone) # Remove negative zone indicator for southern hemisphere
    ko   = 0.9996e0    # UTM scale factor
    lat1 = N/ko/a
    dlat = 1e0
    iterations = 0
    while abs(dlat) > 1e-12 and iterations < MAX_UTM_ITERATIONS:
        A0=1e0-(e2/4e0)-(e22*3e0/64e0)-(e23*5e0/256e0)-(e24*175e0/16384e0)
        A2=(3e0/8e0)*( e2+(e22/4e0)+(e23*15e0/128e0)-(e24*455e0/4096e0) )
        A4=(15e0/256e0)*( e22+(e23*3e0/4e0)-(e24*77e0/128e0) )
        A6=(35e0/3072e0)*( e23-(e24*41e0/32e0) )
        A8=-(315e0/131072e0)*e24
        f1=a*( A0*lat1-A2*math.sin(2e0*lat1)+A4*math.sin(4e0*lat1)-A6*math.sin(6e0*lat1)+A8*math.sin(8e0*lat1) )-N/ko
        f2=a*( A0-2e0*A2*math.cos(2e0*lat1)+4e0*A4*math.cos(4e0*lat1)-6e0*A6*math.cos(6e0*lat1)+8e0*A8*math.cos(8e0*lat1) )
        dlat=-f1/f2
        lat1=lat1+dlat
        iterations += 1
    if iterations >= MAX_UTM_ITERATIONS:
      raise RuntimeError('[ERROR] utm2ell failed to converge after 100 iterations')
    RN  = ell.N(lat1)
    RM  = ell.M(lat1)
    h2  = e2*pow(math.cos(lat1),2e0)/(1e0-e2)
    t   = math.tan(lat1)
    t2  = pow(t,2e0)
    t4  = t2*t2
    t6  = t4*t2
    h22 = pow(h2,2e0)
    h23 = h22*h2
    h24 = h23*h2

    E0  = E/ko/RN
    E1  = E0
    E2  = pow(E0,3e0)/6e0*(1e0+2e0*t2+h2)
    E3  = pow(E0,5e0)/120e0*(5e0+6e0*h2+28e0*t2-3e0*h22+8e0*t2*h2+24e0*t4-4e0*h23+4e0*t2*h22+24e0*t2*h23)
    E4  = pow(E0,7e0)/5040e0*(61e0 + 662e0*t2 + 1320e0*t4 + 720e0*t6)
    lon = (1e0/math.cos(lat1))*(E1-E2+E3-E4)+lcm

    E0 = E/ko
    N1 = (t*pow(E0,2e0))/(2e0*RM*RN)
    N2 = (t*pow(E0,4e0))/(24e0*RM*pow(RN,3e0))*(5e0+3e0*t2+h2-4e0*h22-9e0*h2*t2)
    N3 = (t*pow(E0,6e0))/(720e0*RM*pow(RN,5e0))*(61e0-90e0*t2+46e0*h2+45e0*t4-252e0*t2*h2-5e0*h22+100e0*h23-66e0*t2*h22-90e0*t4*h2+88e0*h24+225e0*t4*h22+84e0*t2*h23-192e0*t2*h24)
    N4 = (t*pow(E0,8e0))/(40320e0*RM*pow(RN,7e0))*(1385e0+3633e0*t2+4095e0*t4+1575e0*t6)
    lat= lat1-N1+N2-N3+N4
    return lat, lon

def ell2utm(lat, lon, ell=Ellipsoid("wgs84"), lcm=None):
    """Ellipsoidal coordinates to UTM.

        Convert ellipsoidal coordinates (actualy longtitude and latitude) to
        UTM coordinates (aka easting and northing).
        If zone is passed in, then it is used for the computation; else, zone
        is computed within the function. The actual zone value used within the
        function is returned in the returned tuple.

        Args:
            lat (float): latitude in radians
            lon (float): longtitude in radians
            ell (Ellipsoid): ellipsoid of choice
            lcm : optional non-standard central meridian in radians

        Returns:
            tuple (float, float, int, int): a tuple of type:
                Northing, Easting, Zone, lcm
    """
    f  = ell.f
    a  = ell.a
    e2 = ell.eccentricity_squared()

    if lcm:
        Zone = 0
    else:
        Zone = floor(math.degrees(lon)/6e0)+31
        Zone = Zone + int(Zone<=0)*60 - int(Zone>60)*60
        lcm = math.radians(Zone*6-183)
    # assert Zone >= 1 and Zone <= 60

    if abs(lat) > math.radians(80):
      print('[WARNING] Latitude outside 80N/S limit for UTM')

    ko = 0.9996e0   # Scale factor
    No = 0e0 if lat > 0e0 else 1e7
    Eo = 500000   # False easting

    lam = lon-lcm
    lam = lam - int(lam>=math.pi) * (math.pi * 2e0)

    RN = ell.N(lat)
    RM = ell.M(lat)

    coslat = math.cos(lat)
    sinlat = math.sin(lat)
    h2     = e2*coslat*coslat/(1e0-e2)
    t      = math.tan(lat)
    n      = f/(2e0-f)

    # powers of various values
    n2   = pow(n,2e0)
    n3   = n2*n
    n4   = n3*n
    t2   = pow(t,2e0)
    t3   = t2*t
    t4   = t3*t
    t6   = t4*t2
    h22  = pow(h2,2)
    h23  = h22*h2
    h24  = h23*h2

    A0 = 1e0+n2/4e0+n4/64e0
    A2 = 3e0/2e0*(n-n3/8e0)
    A4 = 15e0/16e0*(n2-n4/4e0)
    A6 = 35e0/48e0*n3
    A8 = 315e0/512e0*n4
    S  = a/(1e0+n)*(A0*lat-A2*math.sin(2e0*lat)+A4*math.sin(4e0*lat)-A6*math.sin(6e0*lat)+A8*math.sin(8e0*lat))

    E1   = lam*coslat
    E2   = pow(lam,3e0)*pow(coslat,3e0)/6e0*(1e0-t2+h2)
    E3   = pow(lam,5e0)*pow(coslat,5e0)/120e0*(5e0-18e0*t2+t4+14e0*h2-58e0*t2*h2+
               13e0*h22+4e0*h23-64e0*t2*h22-24e0*t2*h23)
    E4   = pow(lam,7e0)*pow(coslat,7e0)/5040e0*(61e0-479e0*t2+179e0*t4-t4*t2)
    E    = Eo+ko*RN*(E1+E2+E3+E4)

    N1 = S/RN
    N2 = pow(lam,2e0)/2e0*sinlat*coslat;
    N3 = pow(lam,4e0)/24e0*sinlat*pow(coslat,3e0)*(5e0-t2+9e0*h2+4e0*h22)
    N4 = pow(lam,6e0)/720e0*sinlat*pow(coslat,5e0)*(61e0-58e0*t2+t4+
            270e0*h2-330e0*t2*h2+445e0*h22+324e0*h23-680e0*t2*h22+
            88e0*h24-600e0*t2*h23-192e0*t2*h24)
    N5 = pow(lam,8e0)/40320e0*sinlat*pow(coslat,7e0)*(1385e0-311e0*t2+543e0*t4-t6)
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
        print('> Testing Point #{} lat={} lon={}'.format(i, degrees(lats[i]), degrees(lons[i])))
        n, e, z, l = ell2utm(lats[i], lons[i], ell)
        print('\tNorthing={} Easting={} Zone={} Central Mer.={}'.format(n, e, z, l))
        print('\tOctave (abs) diffs: dN{} dE{} dZ{} dCM{}'.format(abs(n-octave[i][0]), abs(e-octave[i][1]), abs(z-octave[i][2]), abs(l-octave[i][3])))
        clat, clon = utm2ell(e, n, z, ell)
        if abs(clat-lats[i])>5e-10 or abs(clon-lons[i])>5e-10:
            print('\tERROR Too big discrepancies for station #{}'.format(i))
            print('\tdlat={} dlon={} in decimal degrees'.format(degrees(abs(clat-lats[i])), degrees(abs(clon-lons[i]))))
            print('\tdLat={} dLon={} in seconds'.format(degrees(abs(clat-lats[i]))*3600e0, degrees(abs(clon-lons[i]))*3600e0))
            print('\tInput {}, {} output {}, {}'.format(degrees(lats[i]), degrees(lons[i]), degrees(clat), degrees(clon)))
