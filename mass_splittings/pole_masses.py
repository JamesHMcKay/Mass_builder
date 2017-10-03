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


# get the value of Xi
Mass=np.genfromtxt('models/MSSM/input.txt',usecols=[1])
Xi = Mass[11]

print("Xi = " +  str(Xi))

#### plot pole masses 


plt.figure()

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_iterative_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_iterative_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_iterative_3.txt',usecols=[0,1,2,3,4,5])



x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'-',color='red',label="$Q=2 m_t$",linewidth=1.2) #
plt.plot(x,ymid/x,'-',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
plt.plot(x,yupper1/x,'-',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
plt.plot(x,yupper2/x,'-',color='blue',label="$Q= M$",linewidth=1.2) #
plt.plot(x,yupper3/x,'-',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #

# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='green',alpha=0.5)


if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_explicit_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_explicit_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/pole_mass_n_explicit_3.txt',usecols=[0,1,2,3,4,5])


x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'--',color='red',label="$Q=2 m_t$",linewidth=1.2) #
plt.plot(x,ymid/x,'--',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
plt.plot(x,yupper1/x,'--',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
plt.plot(x,yupper2/x,'--',color='blue',label="$Q= M$",linewidth=1.2) #
plt.plot(x,yupper3/x,'--',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='grey',alpha=0.5)


xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=18)
ylabel(r"$M^0_{\mathrm{pole}} /\hat{M}$",fontsize=18)
plt.xscale('log')

plt.xlim([min(x),5e3])
plt.ylim([0.85,1.10])

plt.tick_params(labelsize=14)

if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_3.pdf")
if (Xi == 1):
	plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_1.pdf")
if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_0.pdf")


plt.figure()

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_iterative_0.txt',usecols=[0,1,2,3,4,5])

if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_iterative_1.txt',usecols=[0,1,2,3,4,5])

if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_iterative_3.txt',usecols=[0,1,2,3,4,5])



x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'-',color='red',label="$Q=2 m_t$",linewidth=1.2) #
plt.plot(x,ymid/x,'-',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
plt.plot(x,yupper1/x,'-',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
plt.plot(x,yupper2/x,'-',color='blue',label="$Q= M$",linewidth=1.2) #
plt.plot(x,yupper3/x,'-',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #


# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='green',alpha=0.5)

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_explicit_0.txt',usecols=[0,1,2,3,4,5])

if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_explicit_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/pole_mass_c_explicit_3.txt',usecols=[0,1,2,3,4,5])



x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]


plt.plot(x,ylower/x,'--',color='red',label="$Q=2 m_t$",linewidth=1.2) #
plt.plot(x,ymid/x,'--',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
plt.plot(x,yupper1/x,'--',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
plt.plot(x,yupper2/x,'--',color='blue',label="$Q= M$",linewidth=1.2) #
plt.plot(x,yupper3/x,'--',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #

# create vector with maximum and minimum values


maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='grey',alpha=0.5)

xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=18)
ylabel(r"$M^+_{\mathrm{pole}} /\hat{M}$",fontsize=18)

plt.xscale('log')
plt.xlim([min(x),5e3])
plt.ylim([0.85,1.10])

plt.tick_params(labelsize=14)

if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_3.pdf")
if (Xi == 1):
	plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_1.pdf")
if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_0.pdf")




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

