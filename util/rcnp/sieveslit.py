import sys
sys.path.append('./util/rcnp')
from sieveslit_raytracer import *
from setfont import *

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,help="Input sieve slit data",default=None)
    args = parser.parse_args()

    if not args.input:
        parser.error("No input file specified.")

    fitter = SieveSlitFit(args.input,x=3,a=1,y=1,b=1)
    fitter.plot_data()
    fitter.global_fit()
    fitter.plot_global_fit()
