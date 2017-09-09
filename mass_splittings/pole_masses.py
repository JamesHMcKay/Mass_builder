#!/usr/bin/python

import numpy as np
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

A=np.genfromtxt('mass_splittings/data/pole_mass_n_iterative.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'-',color='yellow',label="$Q=2 m_t$") #
plt.plot(x,ymid/x,'-',color='red',label="$Q=0.5 m_t$") #
plt.plot(x,yupper1/x,'-',color='green',label="$Q=0.5 M$") #
plt.plot(x,yupper2/x,'-',color='blue',label="$Q= M$") #
plt.plot(x,yupper3/x,'-',color='black',label="$Q=2 M$") #

plt.fill_between(x, ylower/x, yupper1/x,color='grey')
plt.fill_between(x, yupper3/x, ymid/x,color='grey')

A=np.genfromtxt('mass_splittings/data/pole_mass_n_explicit.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'--',color='yellow',label="$Q=2 m_t$") #
plt.plot(x,ymid/x,'--',color='red',label="$Q=0.5 m_t$") #
plt.plot(x,yupper1/x,'--',color='green',label="$Q=0.5 M$") #
plt.plot(x,yupper2/x,'--',color='blue',label="$Q= M$") #
plt.plot(x,yupper3/x,'--',color='black',label="$Q=2 M$") #


plt.fill_between(x, ylower/x, yupper1/x,color='lightgrey')
plt.fill_between(x, yupper3/x, ymid/x,color='lightgrey')


xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=14)
ylabel(r"$M^0_p /\hat{M}$",fontsize=14)
plt.xscale('log')

plt.xlim([min(x),5e3])
plt.ylim([0.85,1.10])

if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_3.pdf")
if (Xi == 1):
	plt.annotate(r'Feynman-\'\,t Hooft gauge $(\xi = 1)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_1.pdf")
if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_n_xi_0.pdf")


plt.figure()

A=np.genfromtxt('mass_splittings/data/pole_mass_c_iterative.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'-',color='yellow',label="$Q=2 m_t$") #
plt.plot(x,ymid/x,'-',color='red',label="$Q=0.5 m_t$") #
plt.plot(x,yupper1/x,'-',color='green',label="$Q=0.5 M$") #
plt.plot(x,yupper2/x,'-',color='blue',label="$Q= M$") #
plt.plot(x,yupper3/x,'-',color='black',label="$Q=2 M$") #

plt.fill_between(x, ylower/x, yupper1/x,color='grey')
plt.fill_between(x, yupper3/x, ymid/x,color='grey')


A=np.genfromtxt('mass_splittings/data/pole_mass_c_explicit.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'--',color='yellow',label="$Q=2 m_t$") #
plt.plot(x,ymid/x,'--',color='red',label="$Q=0.5 m_t$") #
plt.plot(x,yupper1/x,'--',color='green',label="$Q=0.5 M$") #
plt.plot(x,yupper2/x,'--',color='blue',label="$Q= M$") #
plt.plot(x,yupper3/x,'--',color='black',label="$Q=2 M$") #

plt.fill_between(x, ylower/x, yupper1/x,color='lightgrey')
plt.fill_between(x, yupper3/x, ymid/x,color='lightgrey')



xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=14)
ylabel(r"$M^+_p /\hat{M}$",fontsize=14)

plt.xscale('log')
plt.xlim([min(x),5e3])
plt.ylim([0.85,1.10])


if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_3.pdf")
if (Xi == 1):
	plt.annotate(r'Feynman-\'\,t Hooft gauge $(\xi = 1)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_1.pdf")
if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
	plt.savefig("mass_splittings/figures/pole_masses_c_xi_0.pdf")



