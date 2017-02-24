import sys
sys.path.append('./util/rcnp')
from sieveslit_raytracer import *
from setfont import *
import ROOT
import IPython
import colormaps as cmaps
import matplotlib.colors as colors

# load grutinizer libraries
ROOT.gSystem.Load("/projects/ceclub/sullivan/cagragr/e441/lib/libAllGrutinizer")
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,help="Input sieve slit data",default=None)
    args = parser.parse_args()

    if not args.input:
        parser.error("No input file specified.")

    fitter = SieveSlitFit(args.input,x=3,a=1,y=1,b=1)
    #fitter.plot_data()
    fitter.global_fit()



    tfile = ROOT.TFile("./run4022_py.root")
    ttree = tfile.Get("EventTree")


    As = []
    Bs = []
    for i in range(0,ttree.GetEntries()):
        ttree.GetEntry(i)
        rcnp = ttree.TGrandRaiden.GetGRHit(0)
        if rcnp.GR_RAYID(0) != 0:
            continue
        #data.append([rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0),rcnp.GR_PH(0)])
        As.append(fitter.fit_a(rcnp.GR_X(0),rcnp.GR_TH(0)))
        Bs.append(fitter.fit_b(rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0),rcnp.GR_PH(0)))
    bin_contents, x_edges, y_edges = np.histogram2d(As,Bs,bins=500,range=[[-20,20],[-65,65]])

    ## plotting ##
    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    setticks(axes,xlog=True,ylog=True)
    mat = axes.pcolormesh(x_edges,y_edges,bin_contents.T,norm=colors.LogNorm(vmin=max(bin_contents.min(),1),vmax=bin_contents.max()))
    fig.colorbar(mat)
    fitter.plot_global_fit(lazy=True)


    plt.show()
    IPython.embed()
