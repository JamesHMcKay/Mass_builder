#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016
 
#--- plot_MSSM.py ---
 
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



e = 0.3134

pi = 4.0*arctan(1.0)

mw = 80.385
mz = 91.1876

CW = mw/mz
CW2 = power(CW,2)
sw2 = 1.0-CW2
sw = power(sw2,0.5)

c0 = (e**2) * (mw - CW2*mz) / (8. * pi * sw2)

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('models/MSSM/output/mass_splittings.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000





plt.plot(x,y2,'-',color='black',label='Explicit pole mass') #

plt.plot(x,y1,'--',color='black',label='Implicit pole mass') #

plt.plot(x,c0*1000.*ones(size(y1)),'--',color='black',label='Limit') #



xlabel(r"Degenerate mass $M$ (GeV)",fontsize=16)
ylabel(r"$\Delta M$ (Mev)",fontsize=16)
plt.ylim([0,200])
plt.legend(loc=2)


plt.axvspan((9.4e3-0.47e3), (9.4e3+0.47e3), alpha=0.8, color='green')

plt.plot((9.4e3,9.4e3),(0,700),color='red')

plt.savefig("mass_splittings_MSSM.eps")