#!/usr/bin/python

import numpy as np
import pylab
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

#### plot two-loop mass splitting
#
# to generate Figure in paper
# ./mass_builder -g -m itMSSM -i models/MSSM/lists/groupAB.txt -o
# then run plot_M_2loop_iterative(data)
# and plot_M_2loop_explicit(data)
# and use the following settings in the input
#
# ma 1e-1
# mw 80.385
# mz 91.1876
# alpha 0.00781592
#
# convergence tolerance set to 10-6 in ew_spectrum.cpp and 30 iteration limit
#

plt.figure()

A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_2loop.txt',usecols=[0,1,2,3,4])

x=A[:,0]
Q10=A[:,1]*1000
Q100=A[:,2]*1000
Q500=A[:,3]*1000
Q1000=A[:,4]*1000

plt.plot(x,Q10,'--',color='yellow') #
plt.plot(x,Q100,'--',color='red') #
plt.plot(x,Q500,'--',color='black') #
plt.plot(x,Q1000,'--',color='blue') #



# create vector with maximum and minimum values

maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:5])
	minimum[i] = min(A[i,1:5])

plt.fill_between(x, minimum*1000,maximum*1000,color='grey',alpha=0.4)



A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_2loop.txt',usecols=[0,1,2,3,4])
x=A[:,0]
Q10=A[:,1]*1000
Q100=A[:,2]*1000
Q500=A[:,3]*1000
Q1000=A[:,4]*1000


plt.plot(x,Q10,'-',color='yellow',label="$Q=10$") #
plt.plot(x,Q100,'-',color='red',label="$Q=100$") #
plt.plot(x,Q500,'-',color='black',label="$Q=500$") #
plt.plot(x,Q1000,'-',color='blue',label="$Q=1000$") #

# create vector with maximum and minimum values

maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	if (max(A[i,1:5]) > 0 ):
		maximum[i] = max(A[i,1:5])
		minimum[i] = min(A[i,1:5])
		
	

plt.fill_between(x, minimum*1000,maximum*1000,color='green',alpha=0.5)


xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_p^+ - M_p^0$  $($MeV$)$",fontsize=18)


plt.xscale('log')
#plt.yscale('log')

#plt.xlim([min(x),5e3])
plt.ylim([0,180])

plt.tick_params(labelsize=14)

#plt.legend(loc='lower right',fontsize=16,frameon=False)

plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(150,200), xytext=(150,200),fontsize=16)
plt.annotate(r"partial two-loop", xy=(150,190), xytext=(150,190),fontsize=16)

plt.savefig("mass_splittings/figures/2_loop_mass_splitting.pdf")


#######################
### make the legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()


pylab.plot(x,Q10,'-',color='yellow',label="$Q=10$") #
pylab.plot(x,Q100,'-',color='red',label="$Q=100$") #
pylab.plot(x,Q500,'-',color='black',label="$Q=500$") #
pylab.plot(x,Q1000,'-',color='blue',label="$Q=1000$") #



# create a second figure for the legend
figLegend = pylab.figure(figsize = (6.8,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left',ncol=5,frameon=False)

for legobj in leg.legendHandles:
    legobj.set_linewidth(3.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//2loop_legend_1.pdf")


#######################
### make the second legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

# plot something, it doesn't matter what
pylab.plot(x,Q10,'-',color='green',alpha=0.5,label="iterative") #
pylab.plot(x,Q100,'-',color='grey',alpha=0.4,label="non-iterative") #


# create a second figure for the legend
figLegend = pylab.figure(figsize = (6.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False,fontsize=16)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//2loop_legend_2.pdf")



