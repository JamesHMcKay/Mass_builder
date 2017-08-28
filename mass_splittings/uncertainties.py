#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()

A=np.genfromtxt('mass_splittings/uncertainties.txt',usecols=[0,1,2,3,4])
x=A[:,0]
iterative_lower=A[:,1]
iterative_upper=A[:,2]
non_iterative_lower=A[:,3]
non_iterative_upper=A[:,4]

#plt.plot(x,iterative_lower,'--',color='grey')
#plt.plot(x,iterative_upper,'--',color='grey')

#plt.plot(x,non_iterative_lower,'--',color='red')
#plt.plot(x,non_iterative_upper,'--',color='red')

plt.fill_between(x, iterative_lower, iterative_upper,color='grey')

#plt.fill_between(x, non_iterative_lower, non_iterative_upper,color='black')

# plot curves over top of uncertainty bands


A=np.genfromtxt('mass_splittings/mass_splittings_iterative.txt',usecols=[0,1,2,3,4,5])
x2=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x2,ylower,'-',color='yellow',label="$Q=2 m_t$") #
plt.plot(x2,ymid,'-',color='red',label="$Q=0.5 m_t$") #
plt.plot(x2,yupper1,'-',color='green',label="$Q=0.5 M$") #
plt.plot(x2,yupper2,'-',color='blue',label="$Q=1 M$") #
plt.plot(x2,yupper3,'-',color='black',label="$Q=2 M$") #


#A=np.genfromtxt('../Figures/data/mass_splittings_explicit.txt',usecols=[0,1,2,3,4,5])
#x2=A[:,0]
#ylower=A[:,1]
#ymid=A[:,2]
#yupper1=A[:,3]
#yupper2=A[:,4]
#yupper3=A[:,5]


#plt.plot(x2,ylower,'-',color='yellow',label="$Q=2 m_t$") #
#plt.plot(x2,ymid,'-',color='red',label="$Q=0.5 m_t$") #
#plt.plot(x2,yupper1,'-',color='green',label="$Q=0.5 M$") #
#plt.plot(x2,yupper2,'-',color='blue',label="$Q=1 M$") #
#plt.plot(x2,yupper3,'-',color='black',label="$Q=2 M$") #



xlabel(r"Degenerate mass (GeV)",fontsize=14)
ylabel(r"$\Delta M$ (MeV)",fontsize=14)

plt.xscale('log')
plt.xlim([min(x),max(x)])
plt.ylim([0,0.25])

plt.legend(loc='lower left',fontsize=12)

plt.savefig("mass_splittings/uncertainties.eps")
