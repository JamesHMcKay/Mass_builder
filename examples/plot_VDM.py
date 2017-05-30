#!/usr/bin/python

# Mass Builder
 
# James McKay
# May 2017
 
#--- plot_VDM.py ---
 
# plot mass splitting for VDM model


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

A=np.genfromtxt('models/VDM/output/mass_splittings.txt',usecols=[0,1,2,3])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000

plt.plot(x,y2,'-',color='grey',label='Explicit pole mass (fixed $\mu$)') #
plt.plot(x,y3,'-',color='black',label='Explicit pole mass ($\mu = \hat{M}$)') #
plt.plot(x,y1,'--',color='black',label='Iterative pole mass') #


plt.plot(x,(0.210901 )*1000.*ones(size(x)),'-.',color='red') #


xlabel(r"$\hat{M}_V$ (GeV)",fontsize=16)
ylabel(r"$\Delta M=M_V^+-M_V^0$ (MeV)",fontsize=16)
plt.ylim([0,310])
plt.xlim([30,100000])
plt.legend(loc=2,fontsize=12)


plt.savefig("mass_splittings.eps")
