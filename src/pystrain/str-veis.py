#! /usr/bin/env python
'''
    strain.py
    calculates translation, rotation, scale & strain

    input:
        horizontal coordinates for 2 periods, or
        horizontal coordinates & differences
            for a set of sites

    output:
        deformation parateters
'''

import scipy, scipy.linalg
import pyproj
import gmtpy

#
# class Strain
#   deformation ellipse parameters of a set of sites
#   construction: coordinates, coordinate differences
class Strain:
    # constructor
    #   read site ids, coordinates & coordinate differences from 2 scipy arrays
    def _add_set(self, ids, data):
        # import names, coordinates & differences
        self.ids = ids
        self.data = data
        self.residuals = scipy.zeros(scipy.shape(self.data))
        # update deformation & strain parameters
        self.center = scipy.mean(self.data, axis = 0)
        self._update()

    # update deformation parameters
    def _update(self):
        self._update_4()
        self._update_6()
        self._update_strain()

    # 4-parameters deformation
    #   Dx = dx      + y rotation + x scale
    #   Dy =      dy - x rotation + y scale
    def _update_4(self):
        # construct system
        Ax = scipy.zeros((len(self.data), 4))
        Ax[:, 0] = 1.0
        Ax[:, 2] = self.data[:, 0] - self.center[0]
        Ax[:, 3] = self.data[:, 1] - self.center[1]
        Ay = scipy.zeros((len(self.data), 4))
        Ay[:, 1] = 1.0
        Ay[:, 2] = self.data[:, 1] - self.center[1]
        Ay[:, 3] = -self.data[:, 0] + self.center[0]
        A = scipy.concatenate((Ax, Ay), axis = 0)
        del Ax, Ay
        b = scipy.concatenate((self.data[:, 2], self.data[:, 3]))
        # solve for parameters
        parameters, residual, rank, sigma = scipy.linalg.lstsq(A, b)
        self.dx = parameters[0]
        self.dy = parameters[1]
        self.r = parameters[3]
        self.m = parameters[2]
        del parameters
        # compute residuals
        self.residuals[:, 0] = self.data[:, 2] - self.dx - self.m * (self.data[:, 0] - self.center[0]) - self.r * (self.data[:, 1] - self.center[1])
        self.residuals[:, 1] = self.data[:, 3] - self.dy + self.r * (self.data[:, 0] - self.center[0]) - self.m * (self.data[:, 1] - self.center[1])

    # 6-parameters deformation
    #   Dx = tx +      x exx + y exy
    #   Dy =      ty +                 x eyx + y eyy
    def _update_6(self):
        # construct system
        Ax = scipy.zeros((len(self.data), 6))
        Ax[:, 0] = 1.0
        Ax[:, 2] = self.data[:, 0] - self.center[0]
        Ax[:, 3] = self.data[:, 1] - self.center[1]
        Ay = scipy.zeros((len(self.data), 6))
        Ay[:, 1] = 1.0
        Ay[:, 4] = self.data[:, 0] - self.center[0]
        Ay[:, 5] = self.data[:, 1] + self.center[1]
        A = scipy.concatenate((Ax, Ay), axis = 0)
        del Ax, Ay
        b = scipy.concatenate((self.data[:, 2], self.data[:, 3]))
        # solve for parameters
        parameters, residual, rank, sigma = scipy.linalg.lstsq(A, b)
        self.tx = parameters[0]
        self.ty = parameters[1]
        self.exx = parameters[2]
        self.exy = parameters[3]
        self.eyx = parameters[4]
        self.eyy = parameters[5]
        del parameters
        # compute residuals
        self.residuals[:, 2] = self.data[:, 2] - self.tx - self.exx * (self.data[:, 0] - self.center[0]) - self.exy * (self.data[:, 1] - self.center[1])
        self.residuals[:, 3] = self.data[:, 3] - self.ty - self.eyx * (self.data[:, 0] - self.center[0]) - self.eyy * (self.data[:, 1] - self.center[1])

    # update strain parameters
    #   e: total rotation
    #   k: mean scale
    #   deformation ellipse parameters (strain, k_max, k_min, az)
    def _update_strain(self):
        self.e = (self.exy - self.eyx) / 2
        self.k = (self.exx + self.eyy) * 1000000. / 2
        self.strain = scipy.sqrt((self.exx - self.eyy) * (self.exx - self.eyy) + (self.exy + self.eyx) * (self.exy + self.eyx)) * 1000000.
        self.k_max = self.k + self.strain / 2
        self.k_min = self.k - self.strain / 2
        self.az = scipy.degrees(2 * scipy.arctan2(self.exy + self.eyx, self.eyy - self.exx))

#
# function plot_strain
#   plot sites & deformation ellipse
def plot_strain(data_set, datum_init = '4326'):
    # setup plot parameters
    plot = gmtpy.GMT(config = {'PAPER_MEDIA' : 'a4+', 'GRID_CROSS_SIZE_PRIMARY' : '0.3c'})
    # plot coastline & relief
    plot.grdimage('etopo1_greece.nc', J = 'T24/15c', C = 'sea.cpt', P = True)
    plot.pscoast(R = 'etopo1_greece.nc', J = True, D = 'f', G = 'c')
    plot.grdimage('etopo1_greece.nc', J = True, C = 'topo.cpt')
    plot.pscoast(R = True, J = True, B = '2g1', Q = True)
    # create temporary list for 'psvelo'
    temp_list = []
    if datum_init == '4326':    # epsg: WGS84, latlong
        strain_data = [data_set.center[0], data_set.center[1], data_set.k_max, data_set.k_min, data_set.az + 90.]
    else:
        set_projection = pyproj.Proj(init = "epsg:%s" % (datum_init))
        WGS84 = pyproj.Proj(init = 'epsg:4326')
        center = pyproj.transform(set_projection, WGS84, data_set.center[0], data_set.center[1])
        strain_data = [center[0], center[1], data_set.k_max, data_set.k_min, data_set.az + 90.]
    temp_list.append(strain_data)
    print strain_data
    # plot deformation ellipse
    plot.psvelo(in_rows = temp_list, R = True, J = True, S = 'x3c', W = 'thinnest')
    # save postscript
    plot.save('test.eps')


#
# drive
f = open("sing.vel")
n = []
d = []
for l in f:
    sid, x, y, dy, dx = l[:-1].split()
    n.append(sid)
    d.append(list(map(float, (x, y, dx, dy))))
print "data (x, y, de, dn)\n", scipy.array(d)
ds = Strain()
ds._add_set(n, scipy.array(d))
print "strain:", ds.strain, ds.k_max, ds.k_min, ds.az
print "residuals\n", ds.residuals
plot_strain(ds, '2100')

