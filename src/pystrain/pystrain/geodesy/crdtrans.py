from math import sqrt, atan2, atan, pi
from pystrain import ellipsoid

def top2daz(north, east, up):
    """ Compute azimouth, zenith and distance from a topocentric vector.
    """
    distance = math.sqrt(north*north + east*east + up*up)
    a        = math.atan2(east, north) % (math.pi*2) # normalized [0, 2pi]
    zenith   = math.acos(up/distance);
    return distance, a, zenith

def car2top(xi, yi, zi, xj, yj, zj, ell=Ellipsoid("wgs84")):
    """ Transform a vector expressed in cartesian coordinates to the topocentric,
        local system around point i (i.e. North(i), East(i), Up(i)).
    """
    # Cartesian to ellipsoidal for reference point.
    phi_i, lamda_i, h_i = car2ell(xi, yi, zi, ell)

    # Trigonometric numbers.
    cosf = math.cos(phi_i)
    cosl = math.cos(lambda_i)
    sinf = math.sin(phi_i)
    sinl = math.sin(lambda_i)

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
    """ Ellipsoidal to cartesian coordinates.
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

def car2elll(x, y, z, ell=Ellipsoid("wgs84")):
    """ Cartesian to ellipsoidal coordinates. Reference: Fukushima, T., 
        "Transformation from Cartesian to geodetic coordinates
        accelerated by Halley's method", J. Geodesy (2006), 79(12): 689-693
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
