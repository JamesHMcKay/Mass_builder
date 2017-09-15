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


# for xi = 3 we turn threshold corrections off, it doesn't change
# the plot except for some strange bumps at low M
#
# all plots made using
# ./mass_builder -g -m xiMSSM -i models/xiMSSM/lists/ChiChi.txt
#
# input values are
# mw 80.385
# mz 91.1876
#
# different gauges are plotted by chaning Xi in models/MSSM/input.txt


rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

plt.figure()

# get the value of Xi
Mass=np.genfromtxt('models/MSSM/input.txt',usecols=[1])
Xi = Mass[11]

print("Xi = " +  str(Xi))

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/uncertainties_0.txt',usecols=[0,1,2,3,4])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/uncertainties_1.txt',usecols=[0,1,2,3,4])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/uncertainties_3.txt',usecols=[0,1,2,3,4])	


x=A[:,0]
iterative_lower=A[:,1]*1000
iterative_upper=A[:,2]*1000
non_iterative_lower=A[:,3]*1000
non_iterative_upper=A[:,4]*1000


plt.fill_between(x, iterative_lower, iterative_upper,color='#7fbf7f')


plt.fill_between(x, non_iterative_lower, non_iterative_upper,color='grey')#color='#cc4c33')

# plot curves over top of uncertainty bands

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_3.txt',usecols=[0,1,2,3,4,5])


x2=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000

plt.plot(x2,ylower,'-',color='red',label="$Q=2 m_t$",linewidth=1.2) #
plt.plot(x2,ymid,'-',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
plt.plot(x2,yupper1,'-',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
plt.plot(x2,yupper2,'-',color='blue',label="$Q=1 M$",linewidth=1.2) #
plt.plot(x2,yupper3,'-',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_3.txt',usecols=[0,1,2,3,4,5])


x2=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000


xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_p^+ - M_p^0$  $($MeV$)$",fontsize=18)

plt.tick_params(labelsize=14)

plt.xscale('log')
plt.xlim([min(x),max(x)])
plt.ylim([0,250])





if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_3.pdf")

if (Xi == 1):
	plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_1.pdf")

if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_0.pdf")


#######################
### make the first legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

plt.plot(x2,ylower,'-',color='red',label="$Q=2 m_t$") #
plt.plot(x2,ymid,'-',color='black',label="$Q=m_t/2$")
plt.plot(x2,yupper1,'-',color='yellow',label="$Q=M/2$")#
plt.plot(x2,yupper2,'-',color='blue',label="$Q= M$")
plt.plot(x2,yupper3,'-',color='#ff7f00',label="$Q=2 M$")



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
pylab.plot(x2,ylower,'-',color='#7fbf7f',label="iterative uncertainty") #
pylab.plot(x2,ymid,'-',color='grey',label="non-iterative uncertainty") #

# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False,fontsize=16)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_2.pdf")


