#!./bin/grutinizer-python
import sys
sys.path.append('./util/rcnp')
from sieveslit_raytracer import *
from pretty_matplotlib import *
from setfont import *

if __name__=="__main__":
    del sys.argv[0] # grutinizer/root issue
    setfont()


    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,help="Input sieve slit data",default=None)
    args = parser.parse_args()

    if not args.input:
        parser.error("No input file specified.")

    fitter = SieveSlitFit(args.input,x=2,a=2,y=1)
    fitter.global_fit()
    fitter.plot_global_fit()

    plt.show()
