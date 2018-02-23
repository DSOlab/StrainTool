#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from sys  import float_info
from math import floor

class Grid:
    """
        A very simple Grid class to be used within the StrainTensor project.
        A Grid instance has x- and y- axis limits and step sizes (i.e x_min,
        x_max, x_step, y_min, y_max, y_step).
        It is iterable; when iterating, the instance will return the center
        of each cell, starting from the bottom left corner and ending at the
        top right. The iteration is performed row-wise (i.e.
        > [x0+xstep/2, y0+ystep/2]
        > [x0+xstep/2, y0+ystep/2+ystep]
        > [x0+xstep/2, y0+ystep/2+2*ystep]
        > ...
        > [x0+xstep/2, ymax-ystep/2]
        > [x0+xstep/2+xstep, y0+ystep/2]
        > [x0+xstep/2+xstep, y0+ystep/2+ystep]
        > [x0+xstep/2+xstep, y0+ystep/2+2*ystep]

    """
    def __init__(self, x_min, x_max, x_step, y_min, y_max, y_step):
        self.x_min = x_min
        self.x_max = x_max
        self.x_step= x_step
        self.y_min = y_min
        self.y_max = y_max
        self.y_step= y_step
        self.cxi    = 0      # current x-axis tick / index
        self.cyi    = 0      # current y-axis tick / index
        #  Watch out for the float-to-integer conversion! Python's int is basicaly
        #+ a floor() so, 7.99999999.... will became a '7' not an '8'.
        self.xpts   = int((x_max-x_min) / x_step + .49)
        self.ypts   = int((y_max-y_min) / y_step + .49)
        print('[DEGUB] Grid x-axis details: x_min {}, x_max {}, x_step {}, #pts {}, computed end at {}, diff {}'.format(self.x_min, self.x_max, self.x_step, self.xpts, x_min + self.xpts * x_step, abs(x_min + self.xpts * x_step - x_max)))
        assert x_min + self.xpts * x_step <= x_max and abs(x_min + self.xpts * x_step - x_max) < x_step/float(2)
        print('[DEGUB] Grid y-ayis details: y_min {}, y_max {}, y_step {}, #pts {}, computed end at {}, diff {}'.format(self.y_min, self.y_max, self.y_step, self.ypts, y_min + self.ypts * y_step, abs(y_min + self.ypts * y_step - y_max)))
        assert y_min + self.ypts * y_step <= y_max and abs(y_min + self.ypts * y_step - y_max) < y_step/float(2)
        assert self.xpts > 0 and self.ypts > 0

    def __iter__(self):
        return self

    def xidx2xval(self, idx):
        """ Given an index (on x-axis), return the value at the centre of this
            cell.
            The index represents the number of a cell (starting from zero).
        """
        assert idx >= 0 and idx < self.xpts
        return self.x_min + self.x_step/2e0 + self.x_step*idx
    
    def yidx2yval(self, idx):
        """ Given an index (on y-axis), return the value at the centre of this
            cell.
            The index represents the number of a cell (starting from zero).
        """
        assert idx >= 0 and idx < self.ypts
        return self.y_min + self.y_step/2e0 + self.y_step*idx

    def next(self):
        """ Return the (centre of the) next cell (aka x,y coordinate pair).
        """
        if self.cxi >= self.xpts:
            if self.cyi+1 >= self.ypts:
                raise StopIteration
            self.cxi  = 0
            self.cyi += 1
            return self.x_max - self.x_step/2, self.yidx2yval(self.cyi)
        else:
            self.cxi += 1
            return self.xidx2xval(self.cxi-1), self.yidx2yval(self.cyi)

def generate_grid(sta_lst, x_step, y_step):
    """ Given a list of Stations and x- and y-axis step sizes, compute and
        return a Grid instance. The max/min x and y values (of the Grid) are
        extracted from the station coordinates; if needed, they are adjusted
        so that (xmax-xmin) is divisible (without remainder) with xstep.

        Obsviously, sta_lst coordinates (assesed by .lon and .lat) must match
        the x_step and y_step values respectively (same units and reference
        system).

        Parameters:
        -----------
        sta_lst: (list) list of Station instances. Coordinates of the stations are
                 used, using the 'lon' and 'lat' instance variables. Longtitude
                 values are matched to 'x_step' and latitude values are matched
                 to 'ystep'
        x_step:  (float) value of step for the x-axis of the grid.
        y_step:  (float) value of step for the y-axis of the grid.

        Returns:
        --------
        A Grid instance; the min and max values of the Grid are computed from
        the input coordinates (i.e. the sta_lst input list) and then adjusted
        so that the range (xmax-xmin) is divisible with x_step (same goes for
        the y-axis)

        TODO:
        The lines 
        assert divmod(y_max-y_min, y_step)[1] == 0e0 and
        assert divmod(x_max-x_min, x_step)[1] == 0e0
        may throw dues to rounding errors; how can i fix that?
    """
    # Get min/max values of the stations.
    y_min = float_info.max
    y_max = float_info.min
    x_min = float_info.max
    x_max = float_info.min
    for s in sta_lst:
        if   s.lon > x_max:
            x_max = s.lon
        elif s.lon < x_min:
            x_min = s.lon
        if   s.lat > y_max:
            y_max = s.lat
        elif s.lat < y_min:
           y_min = s.lat
    # Adjust max and min to step.
    print("\t[DEBUG] Region: Easting: {:}/{:} Northing: {:}/{:}".format(x_min, x_max, y_min, y_max))
    s      = float((floor((y_max-y_min)/y_step)+1e0)*y_step)
    r      = s-(y_max-y_min)
    y_min -= r/2
    y_max += r/2
    # assert divmod(y_max-y_min, y_step)[1] == 0e0
    s      = float((floor((x_max-x_min)/x_step)+1e0)*x_step)
    r      = s-(x_max-x_min)
    x_min -= r/2
    x_max += r/2
    # assert divmod(x_max-x_min, x_step)[1] == 0e0
    print("\t[DEBUG] Adjusted Easting : from {} to {} with step={} pts={}".format(x_min, x_max, x_step, (x_max-x_min)/x_step))
    print("\t[DEBUG] Adjusted Northing: from {} to {} with step={} pts={}".format(y_min, y_max, y_step, (y_max-y_min)/y_step))
    # return a Grid instance
    return Grid(x_min, x_max, x_step, y_min, y_max, y_step)
