#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

import ROOT

from pretty_style import load_style

def all_bin_edges(taxis):
    output = [taxis.GetBinLowEdge(1)]
    output.extend(taxis.GetBinUpEdge(i)
                  for i in range(1, taxis.GetNbins()+1))
    return output

def plot_root_1d_hist(hist):
    bin_edges = all_bin_edges(hist.GetXaxis())

    bin_contents = [hist.GetBinContent(1)]
    bin_contents.extend(hist.GetBinContent(i)
                        for i in range(1, hist.GetXaxis().GetNbins()+1))

    return [bin_edges, bin_contents]

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

