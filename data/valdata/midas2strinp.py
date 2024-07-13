#!/usr/bin/env python3

import sys
import os
import math
import pandas as pd

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + content)

## add first line
# line_prepender('midas.full', 'code lon lat ve vn sve svn sne tp')

## Read midas
## input data
data = pd.read_csv("midas.full", sep=" ")

data.loc[data['lon'] < -180, 'lon2'] = data['lon'] + 360
data.loc[data['lon'] >= -180, 'lon2'] = data['lon']
#data['lon2'] = data['lon']+360
data['ve2'] = data['ve']*1000.
data['vn2'] = data['vn']*1000.
data['sve2'] = data['sve']*1000.
data['svn2'] = data['svn']*1000.

header=['code', 'lon2', 'lat', 've2', 'vn2', 'sve2', 'svn2', 'sne', 'tp' ]
#data.to_csv('midas.vel', columns=header, index=False, header=False, float_format='%.10f', sep=' ')


#data001 = data.query('(lon2 >= 20) & (lon2 <= 44) & (lat >= 34) & (lat <= 45)')
#data001.to_csv('midas001.vel', columns=header, index=False, header=False, float_format='%.7f', sep=' ')

data002 = data.query('(lon2 > -5 ) & (lon2 <= 5)')
data002 = data002.query('(lat >= 40) & (lat <= 45)')
data002.to_csv('midas002.vel', columns=header, index=False, header=False, float_format='%.10f', sep=' ')

#data003 = data.query('(lon2 >= 94) & (lon2 <= 104) & (lat >= -4) & (lat <= 6)')
#data003.to_csv('midas003.vel', columns=header, index=False, header=False, float_format='%.10f', sep=' ')

#data004 = data.query('(lon2 >= 297) & (lon2 <= 312) & (lat >= -40) & (lat <= -30)')
#data004.to_csv('midas004.vel', columns=header, index=False, header=False, float_format='%.10f', sep=' ')
