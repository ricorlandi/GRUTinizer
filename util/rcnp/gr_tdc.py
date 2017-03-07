#!/usr/bin/python
import sys
sys.path.append('./util/rcnp')
from setfont import *
import ROOT
import IPython
import argparse
import colormaps as cmaps
import matplotlib.colors as colors
import matplotlib.pyplot as plt



# load grutinizer libraries
ROOT.gSystem.Load("/projects/ceclub/sullivan/cagragr/e441/lib/libAllGrutinizer")
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)


def get_program_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,help="Input TDCR histograms (.root) file",default=None)
    parser.add_argument("--output", type=str,help="Output gr_wtdc.hst path",default='./grtdc.hst')
    parser.add_argument("--directory", type=str,help="TDirectory in which histograms are stored",default="")
    parser.add_argument("--fraction_of_mean", type=float,help="Fraction of mean of the uniform part of the timing spectra to use as the lower and upper bounds of the timing spectra",default=0.1)
    args = parser.parse_args()
    if not args.input:
        parser.error("No input file specified.")
    return args;


if __name__=="__main__":
    args = get_program_arguments()
    tfile = ROOT.TFile(args.input)
    hists = [
        tfile.Get(args.directory+"GR_TDCR_X1"),
        tfile.Get(args.directory+"GR_TDCR_U1"),
        tfile.Get(args.directory+"GR_TDCR_X2"),
        tfile.Get(args.directory+"GR_TDCR_U2")
    ]
    # check if histograms were successfully found
    for hist in hists:
        if hist == None:
            print "\nERROR: Histogram not found in root file. If the histogram is \nwithin a TDirectory, you must specify this directory name. E.g. \n\npython util/rcnp/gr_tdc.py hist4041.root --directory /GR_hists\n"
            exit()
    print

    # find boundaries of TDC histogram
    percent = args.fraction_of_mean # use 10% of mean to determine edges of timing spectra
    low,high = (180,200) # range to determine mean from uniform part of tdc spectrum
    for hist in hists:
        mean = 0.
        for i in range(low,high):
            mean += hist.GetBinContent(i)
        mean /= (high-low)
        bound_val = mean*percent

        lower_bound = 0
        upper_bound = 0
        for i in range(1,hist.GetNbinsX()+1):
            counts = hist.GetBinContent(i)

            if lower_bound == 0:
                if counts > bound_val:
                    lower_bound = i-1
            else:
                if upper_bound == 0 and counts < bound_val:
                    upper_bound = i
                    break
        print "6  def   GR_WTDC_"+hist.GetName()[-2:]+" = {"+"{},{}".format(lower_bound,upper_bound)+"}"



    # write out the histogram contents
    output = file(args.output,'w')
    for i in range(1,hists[0].GetNbinsX()+1):
        line = [ str(('%f' % item.GetBinContent(i)).rstrip('0'))  for item in hists]
        line.append("\n")
        output.write("%8s%8s%8s%8s%1s" % tuple(line))
    output.close()
    print "Histogram contents written to", args.output
    print "Replace GRAnalyzer/gr_wtdc.hst with this output file, and add the above printed lines to your hist.def file.\n"
