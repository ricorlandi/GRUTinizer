import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import argparse
# http://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
# http://stackoverflow.com/questions/16333296/how-do-you-create-nested-dict-in-python
import collections
# http://matplotlib.org/examples/pylab_examples/demo_tight_layout.html
import matplotlib.gridspec as gridspec
# https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.optimize.curve_fit.html
from scipy.optimize import curve_fit




def Afit(var,a0,a1,a2):
    x,a = var
    return a0+a1*x+a2*a


class SieveSlitFit(object):
    def __init__(self,filepath,x=1,a=1,y=1):
        self.drawing_holes = []
        self.data_holes = []
        self.A = []
        self.B = []
        self.x = []
        self.exclusions = set()
        self.excluded_holes = []
        self.parse(filepath)
        SieveSlitFit.degree_x = x
        SieveSlitFit.degree_a = a
        SieveSlitFit.degree_y = y

    def parse(self,filepath):
        B = set()
        A = set()
        x = set()

        drawing = False
        data = False
        hole_count = 0
        for line in open(filepath):
            elements = line.split()
            if len(elements) == 0:
                drawing = False
                data = False
                continue
            if '#' in line:
                continue
            if 'Drawing' in line:
                drawing = True
                continue
            if 'Data' in line:
                data = True
                continue

            if '!' in line:
                # mark for exclusion
                row = [float(k) for k in elements[:-1]]
                self.exclusions.add(tuple(row))
                if drawing:
                    self.excluded_holes.append(hole_count)
                    hole_count+=1
                elements = elements[:-1]
            row = [float(k) for k in elements]


            if drawing:
                self.drawing_holes.append(row)
                A.add(row[0])
                B.add(row[1])

            if data:
                self.data_holes.append(row)
                x.add(row[0])
        self.drawing_holes = np.asarray(self.drawing_holes)
        self.data_holes = np.asarray(self.data_holes)
        self.A = sorted(A)
        self.B = sorted(B)
        self.x = sorted(x)

    def plot_y_dependence(self):
        self.fit_params = [[]]*len(self.x)
        for i,x in enumerate(self.x):
            scatter_a = []
            scatter_y = []
            for row in self.data_holes:
                if row[0] == x:
                    scatter_a.append(row[1])
                    scatter_y.append(row[2])
            a_prev = scatter_a[0]
            ys = []
            for j,a in enumerate(scatter_a):
                if a != a_prev or j == len(scatter_a)-1:
                    if j == len(scatter_a)-1:
                        ys.append(scatter_y[j])
                    plt.plot(ys,self.B)
                    params = np.polyfit(ys,self.B,SieveSlitFit.degree_y)
                    self.fit_params[i].append(params)
                    plt.plot(ys,[sum([params[n]*pow(y,len(params)-n-1) for n in range(0,len(params))]) for y in ys],'r--')
                    plt.show()
                    ys = []
                    a_prev = a
                ys.append(scatter_y[j])


            plt.scatter(scatter_a,scatter_y)
            plt.show()

    def fit_y_dependence(self):
        #self.fit_params = []
        #self.fit_params = dict()
        self.fit_params = collections.defaultdict(dict)

        xpt = []
        apt = []
        for i,x in enumerate(self.x):
            scatter_a = []
            scatter_y = []
            for row in self.data_holes:
                if row[0] == x:
                    scatter_a.append(row[1])
                    scatter_y.append(row[2])
            a_prev = scatter_a[0]
            ys = []
            for j,a in enumerate(scatter_a):
                # many y's for a given a, only fit when a changes
                if a != a_prev or j == len(scatter_a)-1:
                    # also fit the last item in the list
                    if j == len(scatter_a)-1:
                        ys.append(scatter_y[j])
                    params = np.polyfit(ys,self.B,SieveSlitFit.degree_y)
                    self.fit_params[x][a_prev] = params
                    ys = []
                    a_prev = a
                ys.append(scatter_y[j])
        apt = sorted(apt)
        pts = set()
        for row in self.data_holes:
            pts.add((row[0],row[1]))
            #xpt.append(row[0])
            #apt.append(row[1])
        data = [list(t) for t in zip(*sorted(pts))]
        params = []
        for i,x in enumerate(data[0]):
            a = data[1][i]
            params.append(self.fit_params[x][a])
        params = np.asarray(params)

        # interpolate data
        interpolant = interpolate.interp2d(data[0],data[1],params[:,0])
        uniquex_interp = np.arange(-600,600,100)
        uniquea_interp = np.arange(-0.1,0.1,0.015)
        uniquea_interp = uniquea_interp[:len(uniquex_interp)]
        domainx = []
        domaina = []
        for x in uniquex_interp:
            for a in uniquea_interp:
                domainx.append(x)
                domaina.append(a)

        rangez = [interpolant(x,domaina[i])[0] for i,x in enumerate(domainx)]


        fig = plt.figure()
        gs1 = gridspec.GridSpec(1, SieveSlitFit.degree_y+1)
        for i in range(0,SieveSlitFit.degree_y+1):
            ax = fig.add_subplot(gs1[i], projection='3d')
            ax.scatter(data[0],data[1],params[:,i],c='red')
            #ax.scatter(domainx,domaina,rangez)
        plt.show()

    def prepare_global_fit_data(self):
        xvals = []
        avals = []
        yvals = []
        Bvals = []
        Avals = []
        for i,row in enumerate(self.data_holes):
            if tuple(row) in self.exclusions:
                continue
            xvals.append(row[0])
            avals.append(row[1])
            yvals.append(row[2])
            Avals.append(self.A[i%5])
            Bvals.append(self.B[i%5])
        return np.asarray(xvals),np.asarray(avals),np.asarray(yvals),np.asarray(Bvals),np.asarray(Avals)

    def global_fit(self,initial_A_params=None):
        xvals,avals,yvals,Bvals,Avals = self.prepare_global_fit_data()
        data = (xvals,avals,yvals)

        prev_x =None
        scatter_a = []
        scatter_y = []
        scatter_B = []
        for i,x in enumerate(xvals):
            if i==0:
                prev_x = x
            if x == prev_x:
                scatter_a.append(avals[i])
                scatter_y.append(yvals[i])
                scatter_B.append(Bvals[i])
            else:
                prev_x = x
                plt.scatter(scatter_a,scatter_y)
                plt.show()
                scatter_a = []
                scatter_y = []
                scatter_B = []
                scatter_a.append(avals[i])
                scatter_y.append(yvals[i])
                scatter_B.append(Bvals[i])


        popt,pcov = curve_fit(bfit,data,Bvals,p0=[1]*(SieveSlitFit.degree_x+1)*(SieveSlitFit.degree_a+1)*(SieveSlitFit.degree_y+1))
        def local_b_fit(x,a,y):
            return bfit([x,a,y],*popt)
        self.fit_b = local_b_fit
        for row in self.data_holes:
            X=row[:3]
            #print X,bfit(X,*popt)


        holes = []
        for i,x in enumerate(xvals):
            holes.append([x,avals[i]])
        prev_hole = None
        counter = 0
        reduced_holes = []
        for i,hole in enumerate(holes):
            if prev_hole == hole:
                counter+=1
            else:
                prev_hole = hole
                if counter >= 4:
                    reduced_holes.append(hole)
                counter = 1
        data = np.asarray(reduced_holes)
        Avals = Avals[:len(data[:,0])]


        plt.scatter(data[:,0],data[:,1])
        plt.show()

        popt2,pcov2 = curve_fit(afit,(data[:,0],data[:,1]),Avals,p0=initial_A_params)
        def local_a_fit(x,a):
            return afit([x,a],*popt2)
        self.fit_a = local_a_fit


        print("A Fit parameters: ")
        #print()
        for par in popt2:
            print '%.17f ' % par
        print("B Fit parameters: ")
        for par in popt:
            print '%.17f ' % par

    def plot_global_B_fit(self):
        for i,x in enumerate(self.x):
            scatter_a = []
            scatter_y = []
            scatter_B = []
            for row in self.data_holes:
                if row[0] == x:
                    scatter_a.append(row[1])
                    scatter_y.append(row[2])
                    scatter_B.append(self.fit_b(row[0],row[1],row[2]))
            plt.scatter(scatter_a,scatter_y,marker='+')
            plt.scatter(scatter_a,scatter_B,c='red')
            #plt.show()

    def plot_global_fit(self):
        scatter_a = []
        scatter_y = []
        scatter_B = []
        scatter_A = []


        xvals,avals,yvals,Bvals,Avals = self.prepare_global_fit_data()
        for i,x in enumerate(xvals):
            a = avals[i]
            y = yvals[i]
            print x,a,y,self.fit_b(x,a,y),self.fit_a(x,a)
            scatter_B.append(self.fit_b(x,a,y))
            scatter_A.append(self.fit_a(x,a))

        true_holes_a = []
        true_holes_b = []
        for b in Bvals:
            for a in Avals:
                true_holes_a.append(a)
                true_holes_b.append(b)
        plt.scatter(true_holes_a,true_holes_b,c='blue')
        plt.scatter(scatter_A,scatter_B,c='red')
        plt.show()


def bfit(var,*c):
    x,a,y = var
    total = 0
    counter = 0
    for i in range(0,SieveSlitFit.degree_x+1):
        sum1 = 0
        for j in range(0,SieveSlitFit.degree_a+1):
            sum2 = 0
            for k in range(0,SieveSlitFit.degree_y+1):
                sum2 += c[counter]*pow(y,k)
                counter+=1
            sum1+=sum2*pow(a,j)
        total+=sum1*pow(x,i)
    return total


def afit(var,*c):
    x,a = var
    sum1 = 0
    counter = 0
    for i in range(0,SieveSlitFit.degree_x+1):
        sum1 += c[i]*pow(x,i)
        counter = i
    for i in range(1,SieveSlitFit.degree_a+1):
        sum1 += c[i+counter]*pow(a,i)
    return sum1
