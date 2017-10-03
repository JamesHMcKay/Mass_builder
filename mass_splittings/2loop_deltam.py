#!/usr/bin/python

import numpy as np
import pylab
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special
from scipy.interpolate import spline

from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

def smooth(A,pts,myorder):	
	x=A[:,0]
	B = zeros([pts,6]);
	xnew = np.linspace(x.min(),x.max(),pts)
	B[:,1:6] = spline(log10(x),A[:,1:6],log10(xnew),order=myorder)
	B[:,0] = xnew
	return B

plt.figure()

A=np.genfromtxt('mass_splittings/data2/deltam_2loop_it.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000
y4=A[:,4]*1000
y5=A[:,5]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
y1_smooth = spline(log10(x),y1,log10(xnew1))
y2_smooth = spline(log10(x),y2,log10(xnew1))
y3_smooth = spline(log10(x),y3,log10(xnew1))
y4_smooth = spline(log10(x),y4,log10(xnew1))
y5_smooth = spline(log10(x),y5,log10(xnew1))


plt.plot(xnew1,y1_smooth,'-',color='red',label="$Q=2 m_t$",linewidth=1.1,zorder=500) #
plt.plot(xnew1,y2_smooth,'-',color='black',label="$Q=m_t/2$",linewidth=1.1,zorder=500)
plt.plot(xnew1,y3_smooth,'-',color='yellow',label="$Q=\hat{M}/2$",linewidth=1.1,zorder=500)#
plt.plot(xnew1,y4_smooth,'-',color='blue',label="$Q= \hat{M}$",linewidth=1.1,zorder=500)
plt.plot(xnew1,y5_smooth,'-',color='#ff7f00',label="$Q=2 \hat{M}$",linewidth=1.1,zorder=500)


# create vector with maximum and minimum values

A[:,4] = A[:,5]
A = smooth(A,5000,3)

x=A[:,0]

# create vector with maximum and minimum values

maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])*1000
	minimum[i] = min(A[i,1:6])*1000
	
plt.fill_between(x, minimum,maximum,color='red',alpha=0.5,zorder=100)	

#plt.fill_between(xnew1, y1_smooth, y3_smooth,color='#7fbf7f')
#plt.fill_between(xnew1, y2_smooth, y5_smooth,color='#7fbf7f')


#plt.plot(x,y1,'-.',color='red',label="$Q=2 m_t$") #
#plt.plot(x,y2,'-.',color='black',label="$Q=m_t/2$")
#plt.plot(x,y3,'-.',color='yellow',label="$Q=\hat{M}/2$")#
#plt.plot(x,y4,'-.',color='blue',label="$Q= \hat{M}$")
#plt.plot(x,y5,'-.',color='#ff7f00',label="$Q=2 \hat{M}$")


# plot the one-loop uncertainty band
A=np.genfromtxt('mass_splittings/data/uncertainties_1.txt',usecols=[0,1,2,3,4])


x=A[:,0]
iterative_lower=A[:,1]*1000
iterative_upper=A[:,2]*1000
non_iterative_lower=A[:,3]*1000
non_iterative_upper=A[:,4]*1000


plt.fill_between(x, iterative_lower, iterative_upper,color='green',alpha=0.5,zorder=10)




xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$  $($MeV$)$",fontsize=18)

plt.xscale('log')
#plt.yscale('log')


plt.ylim([0.0,250])

plt.tick_params(labelsize=14)


#leg = plt.legend(loc='lower left')

#for legobj in leg.legendHandles:
#    legobj.set_linewidth(2.0)

plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(13,30), xytext=(13,30),fontsize=16)
plt.annotate(r"Iterative mass splittings", xy=(13,14), xytext=(13,14),fontsize=16)
plt.savefig("mass_splittings/figures2/deltam_2loop_iterative.pdf")




#######################
### make the legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

# plot something, it doesn't matter what
pylab.plot(x,x,'-',color='#7fbf7f',label="iterative uncertainty") #
pylab.plot(x,x,'-',color='grey',alpha=0.5,label="non-iterative uncertainty") #

# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False,fontsize=16)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_3.pdf")

