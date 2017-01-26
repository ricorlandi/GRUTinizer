import os
import sys
sys.path.append('./util/rcnp')
from sieveslit_raytracer import *
from pretty_matplotlib import * # imports ROOT
from setfont import *
import colormaps as cmaps
import matplotlib_hist_errors


# load grutinizer libraries
ROOT.gSystem.Load("/projects/ceclub/sullivan/cagragr/e441/lib/libAllGrutinizer")

class run_by_run_corrections(object):
    def __init__(self):
        self.histogram_files = []
        self.analysis_callbacks = []
        self.run_end_callback = None
        self.terminal_callback = None

    def add_hist_file(self,filepath):
        self.histogram_files.append(filepath)

    def add_hist_files(self,filelist):
        self.histogram_files.extend(filelist)

    def register_run_end_callback(self, callback):
        self.run_end_callback = callback

    def register_analysis_callback(self, callback):
        self.analysis_callbacks.append(callback)

    def register_terminal_callback(self, callback):
        self.terminal_callback = callback

    def start(self):
        outputs = DefaultOrderedDict(list)
        for filepath in self.histogram_files:
            tfile = ROOT.TFile(filepath,"old")
            if tfile.IsOpen():
                for i,callback in enumerate(self.analysis_callbacks):
                    return_vals = callback(tfile)
                    if return_vals is not None:
                        outputs[filepath].append(return_vals)
                    else:
                        print("Error: Return value of analysis callback is None, for file "+tfile.GetName())
                if self.run_end_callback is not None:
                    self.run_end_callback(return_vals)
                tfile.Close()
            else:
                print("Error: "+filepath+" Not Found. Skipping.")
        if self.terminal_callback is not None:
            self.terminal_callback(outputs)





def angle_drift(tfile):
    hist = tfile.Get("GR_new/B[A]")
    if hist is None:
        return None
    else:
        # ata fit
        xmin,xmax = [-10.,10.]
        hist.GetXaxis().SetRangeUser(xmin, xmax)
        low = hist.GetYaxis().FindBin(-10.)
        high = hist.GetYaxis().FindBin(10.)
        xproj = hist.ProjectionX("xproj",low,high)
        #basic_gaus_fit_root(xproj,plot=False)
        parx,covx = basic_gaus_fit_mpl(xproj,plot=False)
        hist.GetXaxis().SetRangeUser(0,0) # reset range

        # bta fit
        ymin,ymax = [-20.,30.]
        hist.GetYaxis().SetRangeUser(ymin,ymax)
        low = hist.GetXaxis().FindBin(-10.)
        high = hist.GetXaxis().FindBin(10.)
        yproj = hist.ProjectionY("yproj",low,high)
        pary,covy = basic_gaus_fit_mpl(yproj,plot=False)

        return [parx[1],pary[1]]



if __name__=="__main__":

    directory = '/projects/ceclub/sullivan/cagragr/e441/MPIEventLoop/'
    files = os.listdir(directory)
    rootfiles = [directory+x for x in files if 'hist' in x and '.root' in x]

    corrector = run_by_run_corrections()
    #corrector.add_hist_files([rootfiles[0]])
    corrector.add_hist_files(rootfiles)

    def print_run_mean(return_vals):
        print(return_vals)

    def print_all_runs_mean(outputs):
        aoffsets = []
        boffsets = []
        for key,value in outputs.iteritems():
            print(key+": (A/B) mean = ", value)
            aoffsets.append(value[0][0])
            boffsets.append(value[0][1])
        domain = range(0,len(aoffsets))
        plt.plot(domain,aoffsets)
        plt.plot(domain,boffsets)
        plt.show()

    #corrector.register_run_end_callback(print_run_mean)
    corrector.register_terminal_callback(print_all_runs_mean)
    corrector.register_analysis_callback(angle_drift)
    corrector.start()























################################## python projection kept for historical reasons only
#xbin_edges,ybin_edges,bin_contents = plot_root_2d_hist(hist)
# transposed_subset = []
# for j,col in enumerate(bin_contents.T):
#     yedge = ybin_edges[j]
#     if yedge > -10.0 and yedge < 10.0:
#         transposed_subset.append(col)
# transposed_subset = np.asarray(transposed_subset)
# projection = transposed_subset.sum(axis=0) # axis 1 for untransposed array
# fig,axes = plt.subplots()
# axes.better_step(ybin_edges,projection)
# transposed_subset = []
# for i,row in enumerate(bin_contents):
#     xedge = xbin_edges[i]
#     if xedge > -10.0 and xedge < 10.0:
#         transposed_subset.append(row)
# transposed_subset = np.asarray(transposed_subset)
# projection = transposed_subset.sum(axis=0) # axis 1 for untransposed array
# fig,axes = plt.subplots()
# axes.better_step(xbin_edges,projection)
# plt.show()
