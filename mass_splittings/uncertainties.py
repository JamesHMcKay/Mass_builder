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

plt.figure()

A=np.genfromtxt('mass_splittings/data/uncertainties.txt',usecols=[0,1,2,3,4])
x=A[:,0]
iterative_lower=A[:,1]*1000
iterative_upper=A[:,2]*1000
non_iterative_lower=A[:,3]*1000
non_iterative_upper=A[:,4]*1000

#plt.plot(x,iterative_lower,'--',color='grey')
#plt.plot(x,iterative_upper,'--',color='grey')

#plt.plot(x,non_iterative_lower,'--',color='red')
#plt.plot(x,non_iterative_upper,'--',color='red')

plt.fill_between(x, iterative_lower, iterative_upper,color='grey')

plt.fill_between(x, non_iterative_lower, non_iterative_upper,color='lightgrey')

# plot curves over top of uncertainty bands


A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative.txt',usecols=[0,1,2,3,4,5])
x2=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000

plt.plot(x2,ylower,'-',color='yellow',label="$Q=2 m_t$") #
plt.plot(x2,ymid,'-',color='red',label="$Q=0.5 m_t$") #
plt.plot(x2,yupper1,'-',color='green',label="$Q=0.5 M$") #
plt.plot(x2,yupper2,'-',color='blue',label="$Q=1 M$") #
plt.plot(x2,yupper3,'-',color='black',label="$Q=2 M$") #


A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit.txt',usecols=[0,1,2,3,4,5])
x2=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000


#plt.plot(x2,ylower,'--',color='yellow')
#plt.plot(x2,ymid,'--',color='red')
#plt.plot(x2,yupper1,'--',color='green')
#plt.plot(x2,yupper2,'--',color='blue')
#plt.plot(x2,yupper3,'--',color='black')

#plt.fill_between(x2, ylower, yupper1,color='lightgrey')
#plt.fill_between(x2, yupper3, ymid,color='lightgrey')



xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=14)
ylabel(r"$\Delta M = M_p^+ - M_p^0$  $($MeV$)$",fontsize=14)

plt.xscale('log')
plt.xlim([min(x),max(x)])
plt.ylim([0,250])



# get the value of Xi
Mass=np.genfromtxt('models/MSSM/input.txt',usecols=[1])
Xi = Mass[11]

print("Xi = " +  str(Xi))


if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_3.pdf")

if (Xi == 1):
	plt.annotate(r'Feynman-\'\,t Hooft gauge $(\xi = 1)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_1.pdf")

if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_0.pdf")


#######################
### make the first legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

pylab.plot(x2,ylower,'-',color='yellow',label="$Q = 2 m_t$") #
pylab.plot(x2,ymid,'-',color='red',label="$Q = m_t/2$") #
pylab.plot(x2,yupper1,'-',color='green',label="$Q = M/2$") #
pylab.plot(x2,yupper2,'-',color='blue',label="$Q=  M$") #
pylab.plot(x2,yupper3,'-',color='black',label="$Q = 2 M$") #



# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left',ncol=5,frameon=False)

for legobj in leg.legendHandles:
    legobj.set_linewidth(3.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_1.pdf")



#######################
### make the second legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

# plot something, it doesn't matter what
pylab.plot(x2,ylower,'-',color='grey',label="iterative uncertainty") #
pylab.plot(x2,ymid,'-',color='lightgrey',label="non-iterative uncertainty") #


# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_2.pdf")


