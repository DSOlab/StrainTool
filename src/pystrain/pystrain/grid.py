#! /usr/bin/python
#-*- coding: utf-8 -*-

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
        print '[DEGUB] Grid x-axis details: x_min {}, x_max {}, x_step {}, #pts {}, computed end at {}, diff {}'.format(self.x_min, self.x_max, self.x_step, self.xpts, x_min + self.xpts * x_step, abs(x_min + self.xpts * x_step - x_max))
        assert x_min + self.xpts * x_step <= x_max and abs(x_min + self.xpts * x_step - x_max) < x_step/float(2)
        print '[DEGUB] Grid y-ayis details: y_min {}, y_max {}, y_step {}, #pts {}, computed end at {}, diff {}'.format(self.y_min, self.y_max, self.y_step, self.ypts, y_min + self.ypts * y_step, abs(y_min + self.ypts * y_step - y_max))
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
        """ Return the next cell.
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
