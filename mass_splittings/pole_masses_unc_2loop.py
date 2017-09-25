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


A=np.genfromtxt('mass_splittings/data2/pole_mass_unc_1loop_n_Full.txt',usecols=[0,1,2,3,4,5])

A = smooth(A,5000,3)

x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'--',color='red',linewidth=1.2) #
plt.plot(x,ymid/x,'--',color='black',linewidth=1.2) #
plt.plot(x,yupper1/x,'--',color='yellow',linewidth=1.2) #
plt.plot(x,yupper2/x,'--',color='blue',linewidth=1.2) #
plt.plot(x,yupper3/x,'--',color='#ff7f00',linewidth=1.2) #


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='green',alpha=0.5,label="one-loop")

A=np.genfromtxt('mass_splittings/data2/pole_mass_unc_2loop_n_Full.txt',usecols=[0,1,2,3,4,5])

A = smooth(A,5000,3)

x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'-',color='red',linewidth=1.2) #
plt.plot(x,ymid/x,'-',color='black',linewidth=1.2) #
plt.plot(x,yupper1/x,'-',color='yellow',linewidth=1.2) #
plt.plot(x,yupper2/x,'-',color='blue',linewidth=1.2) #
plt.plot(x,yupper3/x,'-',color='#ff7f00',linewidth=1.2) #


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='grey',alpha=0.8,label="two-loop")


A=np.genfromtxt('mass_splittings/data2/pole_mass_unc_2loop_n_AB.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'-.',color='red',linewidth=1.2) #
plt.plot(x,ymid/x,'-.',color='black',linewidth=1.2) #
plt.plot(x,yupper1/x,'-.',color='yellow',linewidth=1.2) #
plt.plot(x,yupper2/x,'-.',color='blue',linewidth=1.2) #
plt.plot(x,yupper3/x,'-.',color='#ff7f00',linewidth=1.2) #


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


#plt.fill_between(x, minimum/x,maximum/x,color='red',alpha=0.2)

plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(110,0.89), xytext=(110,0.89),fontsize=16)


leg = plt.legend(loc='lower left')

for legobj in leg.legendHandles:
    legobj.set_linewidth(0.0)


xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$M^0_p /\hat{M}$",fontsize=18)
plt.xscale('log')

#plt.xlim([min(x),5e3])
plt.ylim([0.85,1.05])

plt.tick_params(labelsize=14)

plt.savefig("mass_splittings/figures2/pole_masses_2loop_Full.pdf")



#####

plt.figure()


A=np.genfromtxt('mass_splittings/data2/deltam_1loop_Full.txt',usecols=[0,1,2,3,4,5])

A = smooth(A,5000,3)

x=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000

# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = 1000*max(A[i,1:6])
	minimum[i] = 1000*min(A[i,1:6])


plt.fill_between(x, minimum,maximum,color='green',alpha=0.7,label='one-loop',zorder=50)



A=np.genfromtxt('mass_splittings/data2/deltam_2loop_AB.txt',usecols=[0,1,2,3,4,5])

A = smooth(A,5000,3)

x=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000

plt.plot(x,ylower,'-.',color='red',linewidth=1.2)
plt.plot(x,ymid,'-.',color='black',linewidth=1.2)
plt.plot(x,yupper1,'-.',color='yellow',linewidth=1.2)
plt.plot(x,yupper2,'-.',color='blue',linewidth=1.2)
plt.plot(x,yupper3,'-.',color='#ff7f00',linewidth=1.2)


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = 1000*max(A[i,1:6])
	minimum[i] = 1000*min(A[i,1:6])


plt.fill_between(x, minimum,maximum,color='red',alpha=0.2,label='partial two-loop')



A=np.genfromtxt('mass_splittings/data2/deltam_2loop_Full.txt',usecols=[0,1,2,3,4,5])

A = smooth(A,5000,3)


x=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = 1000*max(A[i,1:6])
	minimum[i] = 1000*min(A[i,1:6])


plt.fill_between(x, minimum,maximum,color='grey',alpha=0.6,label='two-loop',zorder=50)



plt.xlim([100,1e4])
#plt.ylim([140,177])
plt.ylim([0,250])



xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_p^+ - M_p^0$  $($MeV$)$",fontsize=18)

plt.xscale('log')

plt.tick_params(labelsize=14)

leg = plt.legend(loc='lower right')

for legobj in leg.legendHandles:
    legobj.set_linewidth(0.0)


plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(850,150), xytext=(850,150),fontsize=16)

plt.savefig("mass_splittings/figures2/deltam_2loop_explicit.pdf")



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

