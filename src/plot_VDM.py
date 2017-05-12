#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016
 
#--- plot_example.py ---
 
# Simple routine to plot output from mass splittings example


import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter


fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('models/VDM/output/mass_splittings.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000


plt.plot(x,y2,'-',color='black',label='Explicit pole mass') #

plt.plot(x,y1,'--',color='black',label='Implicit pole mass') #



xlabel(r"MS-bar mass $MV_p$, $MV_0$ (GeV)",fontsize=16)
ylabel(r"$\Delta M=MV_p-MV_0$ (Mev)",fontsize=16)
plt.ylim([-30,30])
#plt.legend(loc=2)


plt.savefig("mass_splittings.eps")

# plot pole masses


fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
#ax.set_yscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('models/VDM/output/masses.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]
y2=A[:,2]


plt.plot(x,y1-x,'-',color='black',label='MVp') #

plt.plot(x,y2-x,'--',color='black',label='MV0') #

#plt.plot(x,x,'--',color='black',label='Tree-level M') #



xlabel(r"MS-bar mass $MV_p$, $MV_0$ (GeV)",fontsize=16)
ylabel(r"Pole masses MV_p, MV_0$ (Gev)",fontsize=16)

plt.legend(loc=2)


plt.savefig("pole_masses.eps")
