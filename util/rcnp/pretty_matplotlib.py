#!/usr/bin/env python

import sys
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from collections import defaultdict, OrderedDict, Callable

import ROOT

from pretty_style import load_style

def all_bin_edges(taxis):
    output = [taxis.GetBinLowEdge(1)]
    output.extend(taxis.GetBinUpEdge(i)
                  for i in range(1, taxis.GetNbins()+1))
    return output

def all_bin_centers(taxis):
    return [taxis.GetBinCenter(i)
                  for i in range(1, taxis.GetNbins()+1)]

def plot_root_1d_hist(hist):
    bin_edges = all_bin_edges(hist.GetXaxis())

    bin_contents = [hist.GetBinContent(i)
                        for i in range(1, hist.GetXaxis().GetNbins()+1)]

    return [bin_edges, bin_contents]

def scatter_from_root_1d_hist(hist):
    bin_centers = all_bin_centers(hist.GetXaxis())

    bin_contents = [hist.GetBinContent(i)
                    for i in range(1, hist.GetXaxis().GetNbins()+1)]

    return [bin_centers, bin_contents]


def plot_root_2d_hist(hist, unfilled_are_blank=True):
    xbins = hist.GetXaxis().GetNbins()
    ybins = hist.GetYaxis().GetNbins()

    xbin_edges = np.array(all_bin_edges(hist.GetXaxis()))
    ybin_edges = np.array(all_bin_edges(hist.GetYaxis()))

    bin_contents = np.empty((xbins,ybins))
    bin_contents[:] = np.nan

    print ybins, xbins
    for i in range(xbins):
        for j in range(ybins):
            content = hist.GetBinContent(i+1, j+1)
            if content or not unfilled_are_blank:
                bin_contents[j,i] = content # Because matplotlib plots things transposed
    bin_contents = np.ma.array(bin_contents, mask=np.isnan(bin_contents))

    return [xbin_edges,
            ybin_edges,
            bin_contents,
        ]

def basic_gaus_fit_root(hist,plot=False):
    if plot:
        res = hist.Fit("gaus","S")
        res.Draw()
        import IPython; IPython.embed()
    else:
        res = hist.Fit("gaus","QSN")
    return np.asarray([res.Value(0),res.Value(1),res.Value(2)]),np.asarray([res.CovMatrix(i,j) for i,j in it.product([0,1,2],[0,1,2])])

def basic_gaus_fit_mpl(hist,plot=False):
    bin_centers,bin_contents = scatter_from_root_1d_hist(hist)
    popt,popv = curve_fit(gaus,bin_centers,bin_contents,p0=[1.,1.,1.])
    if plot:
        fig,axes = plt.subplots()
        bin_edges,bin_contents = plot_root_1d_hist(hist)
        axes.better_step(bin_edges,bin_contents)
        step = (bin_edges[1] - bin_edges[0])/2.0
        domain = np.arange(bin_edges[0],bin_edges[-1]+step,step)
        rng = [gaus(x,*popt) for x in domain]
        plt.plot(domain,rng,color='red')
        plt.show()
    return popt,popv

def gaus(x, *params):
    scale,mu,sigma = params
    return scale/np.sqrt(2*np.pi*sigma)*np.exp(-np.power(x-mu,2)/(2*np.power(sigma,2)))



class DefaultOrderedDict(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))
