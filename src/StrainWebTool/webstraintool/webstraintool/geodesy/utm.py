from math import floor, degrees, radians, pi, sin, cos, tan
import ellipsoid

def ell2utm(lat, lon, ell=Ellipsoid("wgs84"), lcm=None):
    f  = ell.f
    a  = ell.a
    e2 = ell.eccentricity_squared()

    Zone = math.floor(math.degrees(lon)/6)+31
    if Zone <= 0:
        Zone += 60
    else:
        Zone -= 60
    lcm = math.radians(Zone*6-183)

    ko = 0.9996   # Scale factor
    if lat > 0:
        No = 0    # False northing (north)
    else:
        No = 1e7  # False northing (south)
    Eo = 500000   # False easting

    lam = lon-lcm;
    if lam >= math.pi: lam = lam - math.pi*2

    RN = ell.N(lat)
    RM = ell.M(lat)

    coslat = math.cos(lat)
    sinlat = math.sin(lat)
    h2 = e2*coslat*coslat/(1-e2)
    t  = math.tan(lat)
    n  = f/(2-f)

    # powers of various values
    n2   = pow(n,2)
    n3   = n2*n
    n4   = n3*n
    t2   = pow(t,2)
    t3   = t2*t
    t4   = t3*t
    h22  = pow(h2,2)
    h23  = h22*h2
    h24  = h23*h2

    A0 = 1+n2/4+n4/64
    A2 = 3/2*(n-n3/8)
    A4 = 15/16*(n2-n4/4)
    A6 = 35/48*n3
    A8 = 315/512*n4
    S  = a/(1+n)*(A0*lat-A2*math.sin(2*lat)+A4*math.sin(4*lat)-A6*math.sin(6*lat)+A8*math.sin(8*lat))

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
