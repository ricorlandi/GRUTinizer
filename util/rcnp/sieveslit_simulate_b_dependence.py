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
    fitter.global_fit()

    fig = plt.figure()

    first_x = float(fitter.holes[0][0])
    for i in range(0,len(fitter.holes)):
        x = float(fitter.holes[i][0])
        if x != first_x:
            break
        print "Simulating hole "+str(i)
        fitter.simulate_data(hole=i,fig=fig)
    plt.show()


    #IPython.embed()
