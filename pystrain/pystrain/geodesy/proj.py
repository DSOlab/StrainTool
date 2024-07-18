#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from pystrain.geodesy.ellipsoid import Ellipsoid, normalize_angle
import math
import utm
import copy

def ell2utm_list(sta_list_ell):
    # mean longitude and latitude
    mlon = sum([ x.lon for x in sta_list_ell ]) / len(sta_list_ell)
    mlat = sum([ x.lat for x in sta_list_ell ]) / len(sta_list_ell)

    # zone letter and number and letter
    zone_num = utm.latlon_to_zone_number(math.degrees(mlat), math.degrees(mlon))
    zone_let = utm.latitude_to_zone_letter(math.degrees(mlat))

    # perfrom transformations using the same zone number and letter
    sta_list_utm = copy.deepcopy(sta_list_ell)
    for j,s in enumerate(sta_list_utm):
        e, n, zn, zl = utm.from_latlon(math.degrees(s.lat), math.degrees(s.lon), zone_num, zone_let)
        assert(zn==zone_num and zl == zone_let)
        sta_list_utm[j].lon = e
        sta_list_utm[j].lat = n

    return sta_list_utm, zone_num, zone_let

def ell2utm_point(lon, lat, zone_num, zone_let):
    e, n, zn, zl = utm.from_latlon(math.degrees(lat), math.degrees(lon), zone_num, zone_let)
    assert(zn==zone_num and zl == zone_let)
    return e,n
