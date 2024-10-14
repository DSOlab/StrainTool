#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from scipy.interpolate import CloughTocher2DInterpolator
import numpy as np
import copy

def filter_sta_list(sta_list):
    sta_list_gnss = [entry for entry in sta_list if entry.technique.lower() == 'g']
    sta_list_sar = [entry for entry in sta_list if entry.technique.lower() == 's']
    return sta_list_gnss, sta_list_sar

def interpolate_velocity(sta_list_src, sta_list_dest, component):
    if component != 'e' and component != 'n':
        print("[ERROR] Unknown (velocity) component to interpolate! must be either 'n' or 'e'", file=sys.sterr)
        raise RuntimeError
    numpts = len(sta_list_src)
    
    # construct a 2d array with data points
    a = np.zeros(shape=(numpts, 2))
    # construct a 2d array with values
    b = np.zeros(shape=(numpts, 1))

    # populate
    for j, entry in enumerate(sta_list_src):
        a[j,0] = entry.lon
        a[j,1] = entry.lat
        b[j] = entry.ve if component=='e' else entry.vn

    # construct the interpolator instance
    interp = CloughTocher2DInterpolator(a, b)

    # assign and return
    for j in range(len(sta_list_dest)):
        z = interp(sta_list_dest[j].lon, sta_list_dest[j].lat)[0]
        if component=='e': sta_list_dest[j].ve = z
        else: sta_list_dest[j].vn = z

    return sta_list_dest

def interpolate_sar_from_gnss(sta_list_utm):
    sta_list_gnss, sta_list_sar = filter_sta_list(sta_list_utm)
    sta_list_sar = interpolate_velocity(sta_list_gnss, sta_list_sar, 'n')
    return sta_list_gnss + sta_list_sar
