from __future__ import division, print_function
import pdb
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.stats import linregress
from itertools import product as iproduct


class BetterFit:
    def __init__(self, x, y, shift_x=None):
        self.shift_x = 0
        self.shift_y = 0
        self.above = (True, True)
        self.raw_x = None
        self.raw_y = None
        self.x = self.raw_x
        self.y = self.raw_y
        self.set_data(x, y)
        self.set_shift(shift_x)
        self.fit_function = None
        self.minimal_range = 0.15
        self.minimal_points = 20
        self.fits = {}
        self.filtered_keys = []
        self.lines=[]
        self.highlighted_points=[]

    def set_data(self, x, y):
        self.raw_x = np.array(x)
        self.raw_y = np.array(y)
        ind = np.argsort(self.raw_x)
        self.x = self.raw_x[ind]
        self.y = self.raw_y[ind]

    def set_shift(self, shift, axis=0):
        if not shift:
            return
        if axis == 0:
            self.shift_x = shift
        else:
            self.shift_y = shift
        self.shift_data()

    def shift_data(self):
        if self.above[0]:
            self.x -= self.shift_x
        else:
            self.x = self.shift_x - self.x
        if self.above[1]:
            self.y -= self.shift_y
        else:
            self.y = self.shift_y - self.y

    def set_fit_function(self, f):
        self.fit_function = f

    def worth_checking(self, index_range):
        x0, x1 = self.x[np.array(index_range)]
        if np.log10(x1) - np.log10(x0) > self.minimal_range and np.diff(index_range) > self.minimal_points:
            return True
        return False

    def get_fit(self, index_range):
        fit_range = self.x[np.array(index_range)]
        if tuple(index_range) in self.fits:
            return self.fits[tuple(index_range)]
        x = self.x[np.logical_and(self.x > fit_range[0], self.x < fit_range[1])]
        y = self.y[np.logical_and(self.x > fit_range[0], self.x < fit_range[1])]

        try:
            popt, pcov = curve_fit(self.fit_function, x, y)
            yfit = self.fit_function(x, *popt)
            dy = y - yfit
            try:
                err = np.mean(abs(np.log10(y) - np.log10(yfit)))
                if np.isinf(err):
                    raise RuntimeError
            except RuntimeWarning:
                raise RuntimeError
            # slope, intercept, r_value, p_value, std_err = linregress(x,dy)
            deviation = linregress(x, dy)
            this_fit = {
                "success": True,
                "params": popt,
                "params_err": np.sqrt(np.diag(pcov)),
                "error": err,
                "deviation": deviation[0],
                "n_points": len(x),
                "range": fit_range
            }

        except RuntimeError:
            this_fit = {"success": False}
        except TypeError:
            this_fit = {"success": False}

        self.fits[tuple(index_range)] = this_fit
        return this_fit

    def check_ranges(self):
        for ipair in iproduct(range(len(self.x)), repeat=2):
            if self.worth_checking(ipair):
                self.get_fit(ipair)
        self.reset_filters()

    def reset_filters(self):
        self.filtered_keys = [i for i in self.fits if self.fits[i]['success']]

    def filter_fits(self, key, min_val, max_val, getter=None):
        if getter is None:
            getter = lambda x: x
        self.filtered_keys = [i for i in self.filtered_keys if min_val < getter(self.fits[i][key]) < max_val]

    def plot_errors(self):
        import matplotlib.pyplot as plt

        fit_range = [np.diff(np.log10(self.fits[i]['range'])) for i in self.filtered_keys]
        err = [self.fits[i]['error'] for i in self.filtered_keys]
        deviation = [self.fits[i]['deviation'] for i in self.filtered_keys]
        perr = [self.fits[i]['params_err'][1] for i in self.filtered_keys]
        beta = [self.fits[i]['params'][1] for i in self.filtered_keys]

        fig = plt.figure()
        plt.ion()
        ax1 = fig.add_subplot(221)
        line1, = ax1.plot(fit_range, deviation, '.', picker=5)
        ax1.set_title('deviation')
        ax2 = fig.add_subplot(222)
        line2, = ax2.plot(fit_range, beta, '.', picker=5)
        ax2.set_title('beta')
        ax3 = fig.add_subplot(223)
        line3, = plt.plot(fit_range, err, '.', picker=5)
        ax3.set_title('error')
        ax4 = fig.add_subplot(224)
        line4, = ax4.plot(fit_range, perr, '.', picker=5)
        ax4.set_title('params_err')
        self.lines = [line1, line2, line3, line4]
        self.axes = [ax1, ax2, ax3, ax4]
        fig.canvas.mpl_connect('pick_event', self.onpick)
        plt.show()

    def onpick(self, event):
        from matplotlib import pyplot as plt
        if event.artist not in self.lines: return True
        N = len(event.ind)
        if not N: return True
        if self.highlighted_points:
            [i.remove() for i in self.highlighted_points]
            del self.highlighted_points[:]

        figi = plt.figure()
        ax = figi.add_subplot(111)
        ax.loglog(self.x, self.y, '.', alpha=0.6)
        for subplotnum, dataind in enumerate(event.ind):
            for subax in self.axes:
                x=subax.get_lines()[0].get_xdata()[dataind]
                y=subax.get_lines()[0].get_ydata()[dataind]
                scatter_point, = subax.plot(x,y,'or',ms=10)
                self.highlighted_points.append(scatter_point)
            subax.get_figure().canvas.draw()
            self.pprint_fit(self.filtered_keys[dataind])
            fit_range = self.x[np.array(self.filtered_keys[dataind])]
            x = self.x[np.logical_and(self.x > fit_range[0], self.x < fit_range[1])]
            this_fit = self.fits[self.filtered_keys[dataind]]
            ax.loglog(x, self.fit_function(x, *this_fit["params"]))
        figi.show()
        return True

    def pprint_fit(self, fit_index):
        for k, v in self.fits[fit_index].iteritems():
            print("%s\t: " % k)
            print(v)
