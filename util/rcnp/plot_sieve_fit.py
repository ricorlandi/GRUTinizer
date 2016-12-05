#!./bin/grutinizer-python
import sys
sys.path.append('./util/rcnp')
from sieveslit_raytracer import *
from pretty_matplotlib import *
from setfont import *

if __name__=="__main__":
    del sys.argv[0] # grutinizer/root issue
    setfont()


    tfile = ROOT.TFile("hist9026_raytrace.root")
    hist = tfile.Get("Target/B[A]")

    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    setticks(axes,xmajor=10,xminor=2,ymajor=50,yminor=10)
    mat = axes.pcolormesh(*plot_root_2d_hist(hist))
    fig.colorbar(mat)

    # fitting sieveslit
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,help="Input sieve slit data",default=None)
    args = parser.parse_args()

    if not args.input:
        parser.error("No input file specified.")

    fitter = SieveSlitFit(args.input,x=2,a=2,y=1)
    fitter.global_fit()
    fitter.plot_global_fit()




    if False:
        #fit diagnostics
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X = np.arange(-600,600,20)
        A = np.arange(-0.1,0.1,0.0030)[0:len(X)]
        smallX=[]
        smallA=[]
        for i,x in enumerate(X):
            if x>=-450 and x<=350:
                smallX.append(x)
        for i,a in enumerate(A):
            if a>=-0.034 and a<0.064:
                smallA.append(a)
        # print len(smallX),len(smallA)
        X,A = np.meshgrid(smallX,smallA)
        Z = fitter.fit_a(X,A)
        #surf = ax.plot_surface(X,A,Z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0,antialiased=False)
        surf = ax.plot_wireframe(X,A,Z)

        xpts = []
        apts = []
        for row in fitter.data_holes:
            xpts.append(row[0])
            apts.append(row[1])
        ax.scatter(xpts,apts,fitter.fit_a(np.asarray(xpts),np.asarray(apts)),c='red')



    #plt.savefig("/user/sullivan/public_html/sieve_uncor.pdf")
    plt.show()
