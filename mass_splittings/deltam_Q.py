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
from matplotlib import rc

rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')


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

A=np.genfromtxt('mass_splittings/data2/deltam_Q.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000

plt.plot(x,y1,'-',color='green',label='1-loop') #

plt.plot(x,y2,'-',color='red',label='2-loop') #


xlabel(r"renormalisation scale $Q$ (GeV)",fontsize=18)
ylabel(r"$\Delta M$ (Mev)",fontsize=18)
#plt.xlim([50,1000])
#plt.ylim([163,170])
#plt.ylim([149,154])

plt.legend(loc=2)


#plt.axvspan((9.4e3-0.47e3), (9.4e3+0.47e3), alpha=0.8, color='green')
#plt.plot((9.4e3,9.4e3),(0,700),color='red')

plt.savefig("mass_splittings/figures2/mass_splittings_MSSM.eps")
