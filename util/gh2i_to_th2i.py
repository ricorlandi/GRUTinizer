#!/usr/bin/python

import glob
import os
import ROOT
ROOT.gSystem.Load("/projects/ceclub/sullivan/cagragr/e441/lib/libAllGrutinizer")
import sys
sys.path.append('./util/rcnp')
from pretty_matplotlib import *


def convert(input_filename, output_filename):
    tf_in = ROOT.TFile(input_filename)
    tf_out = ROOT.TFile(output_filename, 'RECREATE')
    for key in tf_in.GetListOfKeys():
        obj = key.ReadObj()
        if isinstance(obj, ROOT.GH2I):
            output = ROOT.TH2I(obj)
        elif isinstance(obj, ROOT.GH2D):
            output = ROOT.TH2D(obj)
        elif isinstance(obj, ROOT.GH1D):
            output = explicit_GH1D_to_TH1D(obj)
        else:
            output = obj

        if not isinstance(obj, ROOT.TTree):
            output.Write()

for filename in sys.argv[1:]:
    base = os.path.splitext(filename)[0]
    output = base + '_th2i.root'
    print 'Converting {} to {}'.format(filename, output)
    convert(filename, output)
