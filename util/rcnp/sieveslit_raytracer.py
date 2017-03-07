## sieveslit_raytracer.py

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

import IPython


class SieveSlitFit(object):
    def __init__(self,filepath,x=1,a=1,y=1,b=1):
        self.holes = []
        self.parse(filepath)
        SieveSlitFit.degree_x = x
        SieveSlitFit.degree_a = a
        SieveSlitFit.degree_y = y
        SieveSlitFit.degree_b = b
        self.xvals = None

    def parse(self,filepath):
        for line in open(filepath):
            elements = line.split()
            if len(elements) == 0:
                continue
            if '#' in line:
                continue
            row = [float(k) for k in elements]
            self.holes.append(row)

    def prepare_global_fit_data(self):
        if self.xvals is not None:
            return self.xvals,self.avals,self.yvals,self.bvals,self.Avals,self.Bvals
        xvals = []
        avals = []
        yvals = []
        bvals = []
        Avals = []
        Bvals = []
        for i,row in enumerate(self.holes):
            xvals.append(row[0])
            avals.append(row[1])
            yvals.append(row[2])
            bvals.append(row[3])
            Avals.append(row[4])
            Bvals.append(row[5])
        return np.asarray(xvals),np.asarray(avals),np.asarray(yvals),np.asarray(bvals),np.asarray(Bvals),np.asarray(Avals)

    def plot_data(self):
        xvals,avals,yvals,bvals,Bvals,Avals = self.prepare_global_fit_data()
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


    def global_fit(self,initial_A_params=None):
        if initial_A_params==None: # default parameter value based on static members
            initial_A_params=[1]*((SieveSlitFit.degree_x+1) + (SieveSlitFit.degree_a))

        # perform fit for B
        xvals,avals,yvals,bvals,Bvals,Avals = self.prepare_global_fit_data()
        data = (xvals,avals,yvals,bvals)

        popt,pcov = curve_fit(bfit,data,Bvals,p0=[1]*(SieveSlitFit.degree_x+1)*(SieveSlitFit.degree_a+1)*(SieveSlitFit.degree_y+1)*(SieveSlitFit.degree_b+1))
        def local_b_fit(x,a,y,b):
            return bfit([x,a,y,b],*popt)
        self.fit_b = local_b_fit


        # perform fit for A
        holes = set()
        for i,x in enumerate(xvals):
            holes.add((x,avals[i],Avals[i]))

        data = np.asarray(list(holes))

        popt2,pcov2 = curve_fit(afit,(data[:,0],data[:,1]),data[:,2],p0=initial_A_params)
        def local_a_fit(x,a):
            return afit([x,a],*popt2)
        self.fit_a = local_a_fit


        print("A Fit parameters: ")
        print "Number of params: ",len(popt2)
        for par in popt2:
            print '%.17f ' % par
        print("B Fit parameters: ")
        print "Number of params: ",len(popt)
        for par in popt:
            print '%.17f ' % par

    def plot_global_fit(self,lazy=False):
        scatter_a = []
        scatter_y = []
        scatter_B = []
        scatter_A = []


        xvals,avals,yvals,bvals,Bvals,Avals = self.prepare_global_fit_data()
        for i,x in enumerate(xvals):
            a = avals[i]
            y = yvals[i]
            b = bvals[i]
            scatter_B.append(self.fit_b(x,a,y,b))
            scatter_A.append(self.fit_a(x,a))

        true_holes_a = []
        true_holes_b = []
        for b in Bvals:
            for a in Avals:
                true_holes_a.append(a)
                true_holes_b.append(b)
        plt.scatter(true_holes_a,true_holes_b,c='blue')
        plt.scatter(scatter_A,scatter_B,c='red')
        if not lazy:
            plt.show()

    def plot_reference_holes(self):
        xvals,avals,yvals,bvals,Bvals,Avals = self.prepare_global_fit_data()

        true_holes_a = []
        true_holes_b = []
        for b in Bvals:
            for a in Avals:
                true_holes_a.append(a)
                true_holes_b.append(b)
        plt.scatter(true_holes_a,true_holes_b,c='blue')

    def simulate_data(self,hole=None,fig=None):
        if fig is None:
            fig = plt.figure()

        xvals,avals,yvals,bvals,Bvals,Avals = self.prepare_global_fit_data()


        points = zip(xvals,avals,yvals,bvals)
        data = [xvals,avals,yvals,bvals]
        #covariances = np.cov(data)
        covariances = np.zeros((4,4))
        np.fill_diagonal(covariances,1)

        simulated_data = []

        np.random.seed(99)
        data_holes_a = []
        data_holes_y = []
        if hole is not None:
            for i,x in enumerate(xvals):
                if x != xvals[hole]:
                    continue
                a = avals[i]
                y = yvals[i]
                data_holes_a.append(a)
                data_holes_y.append(y)

        i = 0
        while i < 1000:
            if hole is None:
                index = int(np.random.uniform(0,len(points)))
                simulated_data.append(np.random.multivariate_normal(list(points[index]),covariances))
            else:
                index = hole
                a = avals[index]
                a = np.random.normal(avals[index],0.001)
                y = yvals[index]
                y = np.random.normal(yvals[index],1)
                b = bvals[index]
                b = np.random.normal(bvals[index],0.0005)
                simulated_data.append([xvals[index],a,y,b])
            i += 1

        simulated_data = np.array(simulated_data).T

        scatter_a = []
        scatter_y = []
        scatter_B = []
        scatter_A = []

        axes = fig.add_subplot(2,1,1)

        xvals,avals,yvals,bvals = simulated_data
        for i,x in enumerate(xvals):
            a = avals[i]
            y = yvals[i]
            b = bvals[i]
            scatter_a.append(a)
            scatter_y.append(y)
            scatter_B.append(self.fit_b(x,a,y,b))
            scatter_A.append(self.fit_a(x,a))
        axes.scatter(scatter_a,scatter_y,marker='x',c='red')
        axes.scatter(data_holes_a,data_holes_y,c='blue')


        true_holes_a = []
        true_holes_b = []
        for b in Bvals:
            for a in Avals:
                true_holes_a.append(a)
                true_holes_b.append(b)

        axes = fig.add_subplot(2,1,2)
        axes.scatter(scatter_A,scatter_B,marker='x',c='red')
        axes.scatter(true_holes_a,true_holes_b,c='blue')
        axes.set_ylim(-60,60)
        #axes.set_xlim(-50,50)


    def eval_bfit_given_hole(self,hole_index,x=None,a=None,y=None,b=None):
        xval,aval,yval,bval = [float(x) for x in self.holes[hole_index]]
        if x == None:
            x = xval
        if a == None:
            a = aval
        if y == None:
            y = yval
        if b == None:
            b = bval
        return self.fit_b(x,a,y,b)




def bfit(var,*c):
    x,a,y,b = var
    total = 0
    counter = 0
    for i in range(0,SieveSlitFit.degree_x+1):
        sum1 = 0
        for j in range(0,SieveSlitFit.degree_a+1):
            sum2 = 0
            for k in range(0,SieveSlitFit.degree_y+1):
                sum3 = 0
                for l in range(0,SieveSlitFit.degree_b+1):
                    sum3 += c[counter]*pow(b,l)
                    counter+=1
                sum2+=sum3*pow(y,k)
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


